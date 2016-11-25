#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "cmd_line_opts.h"
#include "build_repeat_families.h"

/* build_repeat_families.c -- build list of repeat families */

// RMH: From Version.c created by the Makefile
extern const char *Version;

int MAXLENGTH;   /* Length of sequence.  It's a misnomer, primarily for backwards compatability. */
#define PADLENGTH 11000 /* should be >= L+MAXOFFSET */
int l;           /* Length of l-mers to look at. Odd or even is OK */
#define smalll 6
#define SEEDMISMATCHES 0 /* allow this many mismatches from seed */
#define HASH_SIZE 16000057 /* prime */
#define SMALLHASH_SIZE 5003 /* prime */
char* SEQUENCE_FILE;      /* Where to get sequence data from. */
char* OUTPUT_FILE;        /* Where to put output. */
char* LMER_TABLE_FILE;    /* Where the lmer table is stored. */
int   MINTHRESH;      /* at end, remove l-mers with freq < MINTHRESH */
int   TANDEMDIST; /* l-mers must be this far apart to avoid tandem repeats */
int   VERBOSE;    /* How chatty?  Real chatty?  Really really chatty?  Super extra really chatty? */

int L;    /* How far to extend. max total length of master is 2*L+l  (L >>> l) */
int MAXOFFSET; /* max offset (5) */
int MAXN;     /* max #occ of lmer (10000) */
int MAXR;     /* max #families (100000) */
int MATCH;     /* match score (1) */
int MISMATCH;  /* mismatch score (-1) */
int GAP;       /* gap (indel) score (-5) */
int CAPPENALTY; /* cap on penalty for exiting alignment (-20) */
int MINIMPROVEMENT; /* amount totalbestscore must improve for each extra letter (3) */
int WHEN_TO_STOP; /* stop if no improvement after extending this far (500) */
#define REMOVE_DEGENERATES 0
#define EXTRAMASK 1 /* default is 1 */
float MAXENTROPY;  /* ignore freq l-mers with entropy greater than this (-0.70) */
int GOODLENGTH; /* minimum length of a good subfamily (50) */

int default_l(int);

/*
 * -------- Internal global variables -------- 
 */
int length; /* length of genome sequence */
char *sequence;  /* [2*MAXLENGTH+3*PADLENGTH]; */
char *removed;   /* [2*MAXLENGTH+3*PADLENGTH]; */

char *master;
char **masters; /* MAXR * (2*L+l) but allocate dynamically */
int *masters_allocated;
int *masterstart; 
int *masterend;
int *pos;
char *rev;  // RMH: Flag indicating if N is a forward/reverse strand lmer
int *upperBoundI;  // RMH: Index into boundaries for lmer upper sequence bound
int N;
int ***score; /* 2 * MAXN * (2*MAXOFFSET+1) */
int **score_of_besty; /* MAXN * (2*MAXOFFSET+1) */
int **maskscore; /* 2 * (2*MAXOFFSET+1) */
  /* Reduce memory use by freeing once no longer needed */
int totalbestscore;       /* alignment score for this y or w */
int besttotalbestscore;   /* best alignment score for any y or w seen so far */
int *bestbestscore;  /* best alignment score for this n for any y or w seen so far */
int *savebestscore;
int besty, bestw, nrepeatocc, nactiverepeatocc, bestnrepeatocc, bestnactiverepeatocc, R;
int bestrepeaty, bestsequencey, bestrepeatw, bestsequencew, repeat2sequence;
int maxrepeaty, minrepeatw;
int **distmat;
int prevbestfreq, prevbesthash;
int *boundaries; // RMH: Fasta sequence boundaries.  This contains
                 //      a non-zero ( zero-based ) position for the
                 //      start of each sequence.  The first sequence
                 //      position is assumed to be 0 and therefore is
                 //      not stored in this array. ( NOTE: This is the
                 //      pre "PADLENGTH" padded position.

#define IUPAC(c) c == 'R' || c=='r' || c=='Y' || c=='y' || c=='M' || c=='m' || c=='K' || c=='k' || c=='W' || c=='w' || c=='S' || c=='s' || c=='B' || c=='b' || c=='D' || c=='d' || c=='H' || c=='h' || c=='V' || c=='v'

void usage() 
{
	fprintf(stderr, "RepeatScout Version %s\n\nUsage: \n"
           "RepeatScout -sequence <seq> -output <out> -freq <freq> -l <l> [opts]\n"
           "     -L # size of region to extend left or right (10000) \n"
           "     -match # reward for a match (+1)  \n"
           "     -mismatch # penalty for a mismatch (-1) \n"
           "     -gap  # penalty for a gap (-5)\n"
           "     -maxgap # maximum number of gaps allowed (5) \n"
           "     -maxoccurrences # cap on the number of sequences to align (10,000) \n"
           "     -maxrepeats # stop work after reporting this number of repeats (10000)\n"
           "     -cappenalty # cap on penalty for exiting alignment of a sequence (-20)\n"
           "     -tandemdist # of bases that must intervene between two l-mers for both to be counted (500)\n"
           "     -minthresh # stop if fewer than this number of l-mers are found in the seeding phase (3)\n"
           "     -minimprovement # amount that a the alignment needs to improve each step to be considered progress (3)\n"
           "     -stopafter # stop the alignment after this number of no-progress columns (100)\n"
           "     -goodlength # minimum required length for a sequence to be reported (50)\n"
           "     -maxentropy # entropy (complexity) threshold for an l-mer to be considered (-.7)\n"
	   "     -v[v[v[v]]] How verbose do you want it to be?  -vvvv is super-verbose.\n",
	   Version );
	exit(1);

}



int main(int argc, char* argv[])
{
  time_t start, finish;
  double duration;
  int x;
  struct llist **headptr;
  struct llist *besttmp;
  FILE *fp;

  start = time(0);
  if( 0 == 1 *
       co_get_string(argc, argv, "-sequence", &SEQUENCE_FILE) *
       co_get_string(argc, argv, "-output", &OUTPUT_FILE) *
       co_get_string(argc, argv, "-freq", &LMER_TABLE_FILE) 
  ) {
    usage();
    exit(1);
  }

  fp = fopen(SEQUENCE_FILE, "ro");
  if( NULL == fp ) {
     fprintf(stderr, "Could not open sequence file %s\n", SEQUENCE_FILE);
     exit(1);
  }
  fseek(fp, 0, SEEK_END);
  MAXLENGTH = ftell(fp);
  fclose(fp);

  co_get_int(argc, argv, "-l", &l) || (l = default_l(MAXLENGTH));

  sequence = (char *) malloc( (2 * MAXLENGTH + 3 * PADLENGTH) * sizeof(char) );
  if( NULL == sequence ) {
    fprintf(stderr, "Could not allocate space for sequence\n");
    exit(1);
  }

  removed = (char *) malloc( (2 * MAXLENGTH + 3 * PADLENGTH) * sizeof(char) );
  if( NULL == removed ) {
    fprintf(stderr, "Could not allocate space for masking array\n");
    exit(1);
  }

  co_get_int(argc, argv, "-L", &L)                           || (L=10000);
  co_get_int(argc, argv, "-match", &MATCH)                   || (MATCH=1);
  co_get_int(argc, argv, "-mismatch", &MISMATCH)             || (MISMATCH=-1);
  co_get_int(argc, argv, "-gap", &GAP)                       || (GAP=-5);
  co_get_int(argc, argv, "-maxgap", &MAXOFFSET)              || (MAXOFFSET=5);
  co_get_int(argc, argv, "-maxoccurrences", &MAXN)           || (MAXN = 10000);
  co_get_int(argc, argv, "-maxrepeats", &MAXR)               || (MAXR = 100000);
  co_get_int(argc, argv, "-cappenalty", &CAPPENALTY)         || (CAPPENALTY=-20);
  co_get_int(argc, argv, "-minimprovement", &MINIMPROVEMENT) || (MINIMPROVEMENT=3);
  co_get_int(argc, argv, "-stopafter", &WHEN_TO_STOP)        || (WHEN_TO_STOP=100);
  co_get_float(argc, argv, "-maxentropy", &MAXENTROPY)       || (MAXENTROPY=-0.7);
  co_get_int(argc, argv, "-goodlength", &GOODLENGTH)         || (GOODLENGTH=50);
  co_get_int(argc, argv, "-tandemdist", &TANDEMDIST)         || (TANDEMDIST=500);
  co_get_int(argc, argv, "-minthresh", &MINTHRESH)           || (MINTHRESH=3);


  VERBOSE = co_get_bool(argc, argv, "-v", &x)    ? 1 :
            co_get_bool(argc, argv, "-vv", &x)   ? 2 :
            co_get_bool(argc, argv, "-vvv", &x)  ? 3 :
            co_get_bool(argc, argv, "-vvvv", &x) ? 20: 0;


  master = (char *)malloc( (2*L+l+1) * sizeof(char));
  master[2*L+l] = (char)NULL;
  if( NULL == master ) {
    fprintf(stderr, "Could not allocate space for master array\n");
    exit(1);
  }

  masters_allocated = (int *)malloc( MAXR * sizeof(int) );
  if( NULL == masters_allocated ) {
    fprintf(stderr, "Could not allocated space for the master's allocation table\n");
    exit(1);
  }

  masterstart = (int *)malloc( MAXR * sizeof(int) );
  if( NULL == masterstart ) {
    fprintf(stderr, "Could not allocated space for the master start positions\n");
    exit(1);
  }

  masterend = (int *)malloc( MAXR * sizeof(int) );
  if( NULL == masterend) {
    fprintf(stderr, "Could not allocated space for the master stop positions\n");
    exit(1);
  }

  // RMH: Initialize the rev array
  rev = (char *)malloc( MAXN * sizeof(char) );
  if( NULL == rev ) {
    fprintf(stderr, "Could not allocated space for the reversed array\n");
    exit(1);
  }

  // RMH: Initialize the upperBoundI array
  upperBoundI = (int *)malloc( MAXN * sizeof(int) );
  if( NULL == upperBoundI ) {
    fprintf(stderr, "Could not allocated space for the upperBoundI array\n");
    exit(1);
  }

  pos = (int *)malloc( MAXN * sizeof(int) );
  if( NULL == pos ) {
    fprintf(stderr, "Could not allocated space for the position array\n");
    exit(1);
  }

  bestbestscore = (int *)malloc( MAXN * sizeof(int) );
  savebestscore = (int *)malloc( MAXN * sizeof(int) );
  if( NULL == bestbestscore || NULL == savebestscore ) {
    fprintf(stderr, "Could not allocated space for internal arrays\n");
    exit(1);
  }

  /* print parameters */
  if(VERBOSE) print_parameters();

  /* build sequence */
  length = build_sequence(sequence,SEQUENCE_FILE);
  // RMH: This is no longer needed.  We use one copy of 
  //      the sequence now.  Saves space and makes further
  //      modifications easier.
  //add_rc(sequence); if(VERBOSE) fprintf(stderr,"Done building sequence\n");
  for(x=0; x<length; x++) removed[x] = 0;

  /* allocate space */
  allocate_space();

  /* build headptr */
  if((headptr = (struct llist **) malloc(HASH_SIZE*sizeof(*headptr))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  if(VERBOSE) fprintf(stderr,"Done allocating headptr\n");
  build_headptr(headptr);
  if(VERBOSE) fprintf(stderr,"Done building headptr\n");

  if( (fp = fopen(OUTPUT_FILE, "w")) == NULL)
  {
    fprintf(stderr,"Could not open input file %s\n", OUTPUT_FILE);  exit(1);
  }

  R = 0;
  prevbestfreq = 1000000000;
  prevbesthash = 0;
  while(1)
  {
    besttmp = find_besttmp(headptr);
    if((besttmp == NULL) || (besttmp->freq < MINTHRESH)) /* 2nd cond added 6/9/04 */
    { 
      if(VERBOSE) printf("Stopped at R=%d since no more frequent l-mers\n",R);
      if(VERBOSE) fprintf(stderr,"Stopped at R=%d since no more frequent l-mers\n",R);
      break;
    }
    /* build master */
    for(x=0; x<l; x++) master[L+x] = sequence[(besttmp->lastocc)+x];
    build_pos(besttmp);      /* computes N */

    // RMH: Added because I was seeing some strange cases.
    if ( N < MINTHRESH )
    {
      printf("Warning: N ( %d ) is less than MINTHRESH ( %d ) yet "
             "besttmp->freq = %d....hmmm\nlmer = ", 
             N, MINTHRESH, besttmp->freq );
      for(x=0; x<l; x++) 
         printf("%c", num_to_char( sequence[(besttmp->lastocc)+x] ) );
      printf("\nh = %d\n", prevbesthash );
      besttmp->freq = N;
      continue;
    }

    if(masters_allocated[R] == 0)
    {
      if((masters[R] = (char *) malloc((2*L+l)*sizeof(*masters[R]))) == NULL)
        { fprintf(stderr,"Out of memory\n");  exit(1); }
      masters_allocated[R] = 1;
    }

    extend_right();
    extend_left();
    if(masterend[R]-masterstart[R]+1 >= GOODLENGTH)
    {
      if(VERBOSE) 
      {
        printf("R=%d: Master lmer hit is ",R);
        for(x=0; x<l; x++) printf("%c",num_to_char(master[L+x]));
        printf("\n");
        fprintf(stderr,"R=%d: Master lmer hit is ",R);
        for(x=0; x<l; x++) fprintf(stderr,"%c",num_to_char(master[L+x]));
        fprintf(stderr,"\n");
      }
      if((SEEDMISMATCHES == 0) && (N != MAXN) && (N != besttmp->freq))
      {
        if(VERBOSE) fprintf(stderr,"WARNING N=%d but goodfreq=%d\n",N,besttmp->freq); /* exit(1); */
        if(VERBOSE) printf("WARNING N=%d but goodfreq=%d\n",N,besttmp->freq); /* exit(1); */
      }
      fprintf(fp,">R=%d\n",R);
      for(x=masterstart[R]; x<=masterend[R]; x++)
      {
        fprintf(fp,"%c",num_to_char(master[x]));
        if((x-masterstart[R])%80 == 79) fprintf(fp,"\n");
      }
      if((x-masterstart[R])%80 > 0) fprintf(fp,"\n");
      if(VERBOSE) 
      {
        printf("R=%d: N is %d\n",R,N);
        fprintf(stderr,"R=%d: N is %d\n",R,N);
        printf("AFTER EXTENDING LEFT, FAMILY %d IS:\n",R);
        printf("%d letters total:\n",masterend[R]-masterstart[R]+1);
        for(x=masterstart[R]; x<=masterend[R]; x++) 
        {
          printf("%c",num_to_char(master[x]));
          if((x-masterstart[R])%80 == 79) printf("\n");
	}
        if((x-masterstart[R])%80 > 0) printf("\n");
        fprintf(stderr,"AFTER EXTENDING LEFT, FAMILY %d IS:\n",R);
        fprintf(stderr,"%d letters total:\n",masterend[R]-masterstart[R]+1);
        for(x=masterstart[R]; x<=masterend[R]; x++) 
	{
          fprintf(stderr,"%c",num_to_char(master[x]));
          if((x-masterstart[R])%80 == 79) fprintf(stderr,"\n");
	}
        if((x-masterstart[R])%80 > 0) fprintf(stderr,"\n");
      }

      finish = time(0);
      duration = difftime(finish,start);
      if(VERBOSE) printf("Time1 to this point is %.1f sec = %.1f min = %.1f hr\n",duration, duration/60.0, duration/3600.0);
      if(VERBOSE) fprintf(stderr,"Time1 to this point is %.1f sec = %.1f min = %.1f hr\n",duration, duration/60.0, duration/3600.0);

      R++;
      if(R == MAXR) break;

      mask_headptr(headptr); /* removed[x]=1 means l-mer x was removed */

      finish = time(0);
      duration = difftime(finish,start);
      if(VERBOSE) printf("Time2 to this point is %.1f sec = %.1f min = %.1f hr\n\n",duration, duration/60.0, duration/3600.0);
      if(VERBOSE) fprintf(stderr,"Time2 to this point is %.1f sec = %.1f min = %.1f hr\n\n",duration, duration/60.0, duration/3600.0);
    }
    else
    {
      R++;
      if(R == MAXR) break;
      mask_headptr(headptr);
      R--;
    }
  }

  finish = time(0);
  duration = difftime(finish,start);
  if(VERBOSE) printf("Program duration is %.1f sec = %.1f min = %.1f hr\n",duration, duration/60.0, duration/3600.0);
  fprintf(stderr,"Program duration is %.1f sec = %.1f min = %.1f hr\n",duration, duration/60.0, duration/3600.0);

  return 0;
}

/* ************************************************************** */

void build_headptr(struct llist **headptr)
{
  FILE *fp;
  char string[l];
  int thisfreq, thisocc, x, h, n;
  struct llist *tmp;

  if( (fp = fopen(LMER_TABLE_FILE, "r")) == NULL)
  {
    fprintf(stderr,"Could not open input file %s\n", LMER_TABLE_FILE);  exit(1);
  }

  /* set each head pointer *headptr to NULL */
  for(h=0; h<HASH_SIZE; h++)
    headptr[h] = NULL;

  n = 0;
  while(1)
  {
    fscanf(fp,"%s",string);
    fscanf(fp,"%d",&thisfreq);
    fscanf(fp,"%d",&thisocc);

    /* printf("Processing string %s thisfreq %d\n",string,thisfreq); */

    if(string[0] < 10) break; /* end of file */
    for(x=0; x<l; x++) string[x] = char_to_num(string[x]); /* must translate */

    h = hash_function(string);
    if(h<0) continue;
    tmp = headptr[h];
    while(tmp != NULL)
    {
      if(lmermatcheither(sequence+(tmp->lastocc),string) == 1) /* forward or rc match */
      {      
        /* already in there.  Means we have reached end of file: all done */
        break;
      }
      tmp = tmp->next;
    }
    /* add new guys */
    if(tmp != NULL) break; /* all done */
    if(compute_entropy(sequence+thisocc) > MAXENTROPY) continue; /* skip */
    if((tmp = (struct llist *) malloc(sizeof(*tmp))) == NULL)
    {
      fprintf(stderr,"Out of memory\n");  exit(1);
    }
    n++;
    tmp->freq = thisfreq;
    tmp->lastocc = thisocc;
    tmp->next = headptr[h]; /* either NULL, or some other l-mer with hash h */
    tmp->pos = NULL;
    headptr[h] = tmp; 
  }

  fclose(fp);
  trim_headptr(headptr);
  build_all_pos(headptr);
}

void trim_headptr(struct llist **headptr)
{
  int h;
  struct llist *tmp, *prevtmp, *nexttmp;

  /* remove l-mers with freq < MINTHRESH */
  for(h=0; h<HASH_SIZE; h++)
  {
    prevtmp = NULL;
    tmp = headptr[h];
    while(tmp != NULL)
    {
      if(tmp->freq >= MINTHRESH)
      { /* don't remove tmp */
        // RMH: Let build_all_pos rebuild this ( based on seq boundaries )
        tmp->freq = 0;
        prevtmp = tmp;
        tmp = tmp->next;
        continue;
      }
      /* remove tmp */
      nexttmp = tmp->next;
      free(tmp);
      tmp = nexttmp;
      if(prevtmp == NULL) /* first guy in linked list */
        headptr[h] = tmp;
      else
        prevtmp->next = tmp;
    }
  }
}

void build_all_pos(struct llist **headptr)
{
  int x, h, pos1, pos2;
  // RMH: Added for rev/fwd strand mods
  struct posllist *postmp2;
  struct llist *tmp;
  struct posllist *postmp, *prevpostmp, *nextpostmp;
  int currBoundary = 0;

  if(VERBOSE) fprintf(stderr,"Starting build_all_pos 1\n");

  for(x=l-1; x<length-l+1; x++)
  {
    // RMH: Make sure we are not on a sequence boundary.  
    //      NOTE: I have decided to adjust the frequency
    //            also at this stage.  This should take care
    //            of build_lmer_table including elements
    //            which cross sequence boundaries.  I don't
    //            like adjusting it here but until build_lmer_table
    //            is revamped or we switch to elmer....
    //
    if ( boundaries[currBoundary] > 0 )
    {
      if ( x == boundaries[currBoundary] + PADLENGTH )
      {
        currBoundary++; 
      }
      if ( x+l > boundaries[currBoundary] + PADLENGTH )
      {
        if(VERBOSE) printf("Skipping lmer @ %d due to sequence boundary\n", (x+l) );
        continue;
      }
    }

    /* look for this l-mer in headptr */
    // RMH: This may be a good place for an optimisation.
    //      We could use the bit table trick to disqualify
    //      lmers that can't be in the hash table. Thus
    //      reducing the hash table lookup cost.
    h = hash_function(sequence+x);
    if(h<0) continue;
    tmp = headptr[h];
    while(tmp != NULL)
    {
      // RMH: Changed this to "either" now that we are not loading
      //      into memory the reverse complemented sequence.
      if( ((SEEDMISMATCHES==0) && (lmermatcheither(sequence+(tmp->lastocc),sequence+x))) /* seq incl rc */ ||
          ((SEEDMISMATCHES>0) && (mismatches(sequence+(tmp->lastocc),sequence+x)<=SEEDMISMATCHES)))
      {
        /* a hit.  Add position x to tmp->pos */
        if((postmp = (struct posllist *) malloc(sizeof(*postmp))) == NULL)
        {
          fprintf(stderr,"Out of memory\n");  exit(1);
        }
        postmp->this = x;
        postmp->next = tmp->pos;
        tmp->pos = postmp;
        tmp->freq++;
      }

      tmp = tmp->next;
    }
  }

  if(VERBOSE) fprintf(stderr,"Starting build_all_pos 2\n");

  /* extra code to remove position pairs within TANDEMDIST */
  // RMH: NOTE - Tandem repeat removal creates assymetry in the
  //             algorithm.  When two lmers are within TANDEMDIST
  //             apart the first one is removed and the last one
  //             kept.  Obviously this is dependent on the strand
  //             you feed the program.  Setting tandemdist to 0 
  //             will restore symmetry for testing purposes.
  for(h=0; h<HASH_SIZE; h++)
  {
    tmp = headptr[h];
    while(tmp != NULL)
    {
      /* mark within-TANDEMDIST guys as negative */
      postmp = tmp->pos;
      while((postmp != NULL) && (postmp->next != NULL))
      {
        pos1 = postmp->this;
        // RMH: Fwd/Rev strand mods
        // Tandem distance needs to be measured from the last
        // occurance of a same-stranded lmer.
        postmp2 = postmp;
        while ( ( postmp2 = postmp2->next ) != NULL )
        {
          pos2 = postmp2->this;
          if ( pos1-pos2 >= TANDEMDIST ) 
            break;
          if ( lmermatch(sequence+pos1,sequence+pos2) )
          { 
            postmp->this = -pos1; 
            break; 
          }
        }
        postmp = postmp->next;
      }

      /* remove guys marked as negative */
      postmp = tmp->pos;
      prevpostmp = NULL;
      while(postmp != NULL)
      {
        if(postmp->this >= 0)
        { /* don't remove postmp */
          prevpostmp = postmp;
          postmp = postmp->next;
          continue;
        }
        /* remove postmp */
        nextpostmp = postmp->next;
        // RMH: Reduce frequency as we remove them?
        //      TODO: Consider if this is what we want.
        tmp->freq--;
        removed[-postmp->this] = 1;
        //
        free(postmp);
        postmp = nextpostmp;
        if(prevpostmp == NULL) /* first guy in linked list */
          tmp->pos = postmp;
        else
          prevpostmp->next = postmp;
      }

      tmp = tmp->next;
    }
  }
  return;

  if(VERBOSE) fprintf(stderr,"Done with  build_all_pos\n");
}

struct llist *find_besttmp(struct llist **headptr)
{
  int h;
  struct llist *tmp, *besttmp;
  int bestfreq;

  /* first, try to match prevbestfreq */
  for(h=prevbesthash; h<HASH_SIZE; h++)
  {
    tmp = headptr[h];
    while(tmp != NULL)
    {
      if(tmp->freq == prevbestfreq)  
      {
        prevbesthash = h;
        return tmp;
      }
      tmp = tmp->next;
    }
  }
  /* otherwise, just find best */
  besttmp = NULL;
  bestfreq = 0;
  for(h=0; h<HASH_SIZE; h++)
  {
    tmp = headptr[h];
    while(tmp != NULL)
    {
      if(tmp->freq > bestfreq)
      {
        besttmp = tmp;
        bestfreq = tmp->freq;
        prevbesthash = h;
      }
      tmp = tmp->next;
    }
  }
  prevbestfreq = bestfreq;

  return besttmp;
}

void allocate_space()
{
  int r, x, n;

  /* char **masters; MAXR * (2*L+l) but allocate dynamically */
  if((masters = (char **) malloc(MAXR*sizeof(*masters))) == NULL)
  { fprintf(stderr,"Out of memory\n");  exit(1); }
  for(r=0; r<MAXR; r++) masters_allocated[r] = 0;

  if((distmat = (int **) malloc((l+1)*sizeof(*distmat))) == NULL)
  { fprintf(stderr,"Out of memory\n");  exit(1); }
  for(x=0; x<=l; x++)
  {
    if((distmat[x] = (int *) malloc((l+1)*sizeof(*distmat[x]))) == NULL)
    { fprintf(stderr,"Out of memory\n");  exit(1); }
  }

  /* int ***score; 2 * MAXN * (2*MAXOFFSET+1) */
  if((score = (int ***) malloc(2*sizeof(*score))) == NULL)
  { fprintf(stderr,"Out of memory\n");  exit(1); }
  for(x=0; x<2; x++)
  {
    if((score[x] = (int **) malloc(MAXN*sizeof(*score[x]))) == NULL)
    { fprintf(stderr,"Out of memory\n");  exit(1); }
    for(n=0; n<MAXN; n++)
    {
      if((score[x][n] = (int *) malloc((2*MAXOFFSET+1)*sizeof(*score[x][n]))) == NULL)
      { fprintf(stderr,"Out of memory\n");  exit(1); }
    }
  }

  /* int **score_of_besty; MAXN * (2*MAXOFFSET+1) */
  if((score_of_besty = (int **) malloc(MAXN*sizeof(*score_of_besty))) == NULL)
  { fprintf(stderr,"Out of memory\n");  exit(1); }
  for(n=0; n<MAXN; n++)
  {
    if((score_of_besty[n] = (int *) malloc((2*MAXOFFSET+1)*sizeof(*score_of_besty[n]))) == NULL)
    { fprintf(stderr,"Out of memory\n");  exit(1); }
  }

  /* int **maskscore; 2 * (2*MAXOFFSET+1) */
  if((maskscore = (int **) malloc(2*sizeof(*maskscore))) == NULL)
  { fprintf(stderr,"Out of memory\n");  exit(1); }
  for(x=0; x<2; x++)
  {
    if((maskscore[x] = (int *) malloc((2*MAXOFFSET+1)*sizeof(*maskscore[x]))) == NULL)
    { fprintf(stderr,"Out of memory\n");  exit(1); }
  }
}

void print_parameters()
{
  printf("Parameters:\n");
  printf("SEQUENCE_FILE %s\n",SEQUENCE_FILE);
  printf("LMER_TABLE_FILE %s\n",LMER_TABLE_FILE);
  printf("MAXLENGTH %d\n",MAXLENGTH);
  printf("PADLENGTH %d\n",PADLENGTH);
  printf("l %d\n",l);
  printf("SEEDMISMATCHES %d\n",SEEDMISMATCHES);
  printf("HASH_SIZE %d\n",HASH_SIZE);
  printf("L %d\n",L);
  printf("MAXOFFSET %d\n",MAXOFFSET);
  printf("MAXN %d\n",MAXN);
  printf("MAXR %d\n",MAXR);
  printf("MATCH %d MISMATCH %d GAP %d\n",MATCH,MISMATCH,GAP);
  printf("MINTHRESH %d\n",MINTHRESH);
  printf("CAPPENALTY %d MINIMPROVEMENT %d\n",CAPPENALTY,MINIMPROVEMENT);
  printf("WHEN_TO_STOP %d\n",WHEN_TO_STOP);
  printf("REMOVE_DEGENERATES %d\n",REMOVE_DEGENERATES);
  printf("TANDEMDIST %d\n",TANDEMDIST);
  printf("EXTRAMASK %d\n",EXTRAMASK);
  printf("MAXENTROPY %f\n",MAXENTROPY);
}

int build_sequence(char *sequence, char *filename)
{
  int i, j, seq;
  char c;
  FILE *fp;
  int boundariesSize = 100;

  // RMH: Initialize the boundaries array
  boundaries = (int *)malloc( boundariesSize * sizeof(int) );
  if( NULL == boundaries ) {
    fprintf(stderr, "Could not allocated space for the boundaries array\n");
    exit(1);
  }
  
  if( (fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr,"Could not open input file %s\n", filename);  exit(1);
  }
  for(i=0; i<PADLENGTH; i++) sequence[i] = 99;
  i = 0;
  seq = 0;
  while( !feof(fp) && (i<MAXLENGTH) )
  { /* process one line of file */
    c = getc(fp);
    if(c == EOF) continue;
    if(c == '\n') continue;
    if(c == '>')
    {
      for(j=0; (getc(fp) != '\n') && !feof(fp); j++)
        ;
      // RMH: Store boundary info
      if ( seq > 0 )
      {
        if ( seq > boundariesSize )
        {
          boundariesSize += 100;
          boundaries = realloc( boundaries, boundariesSize * sizeof(int) );
          if( NULL == boundaries ) {
            fprintf(stderr, "Could not allocated more space for the "
                    "boundaries array\n");
            exit(1);
          }
        }
        //printf("Adding sequence boundary[%d] = %d, %d\n", 
        //                      ( seq-1 ), i, i+PADLENGTH );
        boundaries[seq-1] = i;
      }
      seq++;
    }
    else
    {
      if(c > 64) 
      {
        sequence[i+PADLENGTH] = char_to_num(c);
        i++;
      }
      for(j=0; ((c = getc(fp)) != '\n') && !feof(fp) && (i<MAXLENGTH); j++)
      {
        if(c > 64) 
        {
          sequence[i+PADLENGTH] = char_to_num(c);
          i++;
        }
      }
    }
  }
  if ( seq == 0 )
    seq = 1;
  if ( seq+1 >= boundariesSize )
  {
    boundariesSize += 100;
    boundaries = realloc( boundaries, boundariesSize * sizeof(int) );
    if( NULL == boundaries ) {
      fprintf(stderr, "Could not allocated more space for the "
              "boundaries array\n");
      exit(1);
    }
  }
  boundaries[seq-1] = i + 1;
  //printf("Adding sequence boundary[%d] = %d, %d\n", seq-1, (i + 1), (i+1)+PADLENGTH );
  boundaries[seq] = 0;
  fclose(fp);
  return (i+PADLENGTH); /* length of genome */
}

void add_rc(char *sequence)
{
  int x;

  for(x=0; x<PADLENGTH; x++)
    sequence[length+x] = 99;
  for(x=0; x<length; x++)
  {
    if(sequence[length-1-x] == 99)
      sequence[length+PADLENGTH+x] = 99;
    else
      sequence[length+PADLENGTH+x] = 3 - sequence[length-1-x];
  }
  length = 2*length + PADLENGTH;
}

void build_pos(struct llist *besttmp)
{
  struct posllist *postmp;
  int x = 0;

  postmp = besttmp->pos;
  N = 0;
  while(postmp != NULL)
  {
    if(removed[postmp->this] == 0)
    {
      pos[N] = postmp->this;
      // RMH: Fwd/Rev Mods
      // NOTE: This could be more efficiently calculated in build_all_pos
      if ( lmermatch(sequence+besttmp->lastocc,sequence+postmp->this) )
        rev[N] = 0;
      else 
        rev[N] = 1;
      // RMH: Find boundaries for lmer
      upperBoundI[N] = -1;
      x = 0;
      while ( boundaries[x] != 0 )
      {
        if ( boundaries[x] + PADLENGTH > pos[N] )
        {
          upperBoundI[N] = x;
          break; 
        }
        x++;
      }
      N++;
      if(N == MAXN) break;    
    }
    postmp = postmp->next;
  }
  return;
}

void extend_right()
{
  int y, n, bestscore, tempscore, offset; /* bestscore is alignment score for single n */
  char a, besta;
  int newtotalbestscore, newtotalbestscore_a; /* for best a and that a, resp. */
  besta = 0;

  /* initialize score[(L+l-1)%2][][] */
  //
  //   Initialize running score by initializing the score for
  //   the last base in lmer.  This is simply the seed score with
  //   -MAXOFFSET to +MAXOFFSET gaps.
  //
  y = L+l-1;
  for(n=0; n<N; n++)
  {
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      score[y%2][n][offset+MAXOFFSET] = l * MATCH; 
      if(offset < 0) score[y%2][n][offset+MAXOFFSET] += -offset * GAP;
      if(offset > 0) score[y%2][n][offset+MAXOFFSET] += offset * GAP;
    }
    bestbestscore[n] = l * MATCH;
  }

  /* extend right */
  besty = L+l-1;
  for(n=0; n<N; n++)
  {
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
      score_of_besty[n][offset+MAXOFFSET] = score[y%2][n][offset+MAXOFFSET];
  }
  besttotalbestscore = 0;
  compute_totalbestscore_right(y); /* should be N*l*MATCH */
  if(VERBOSE >= 4) print_totalbestscore_right(y); 
  for(y=L+l; y<2*L+l; y++) /* determine master[y] and score[y%2][][] */
  {
    newtotalbestscore = 0;
    for(a=0; a<4; a++) 
    {
      newtotalbestscore_a = 0;
      for(n=0; n<N; n++)
      {
        bestscore = bestbestscore[n] + CAPPENALTY;
        if(bestscore < 0) bestscore = 0;
        for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
        {
          tempscore = compute_score_right(y,n,offset,a);
          if(tempscore > bestscore)
            bestscore = tempscore;
        }
        newtotalbestscore_a += bestscore;
      }
      if(newtotalbestscore_a > newtotalbestscore)
      {
        newtotalbestscore = newtotalbestscore_a;
        besta = a;
      }
    }
    master[y] = besta;
    for(n=0; n<N; n++)
    {
      for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
        score[y%2][n][offset+MAXOFFSET] = compute_score_right(y,n,offset,besta);
    }
    compute_totalbestscore_right(y); 
    if(VERBOSE >= 4) print_totalbestscore_right(y); 
  
    if(y - besty >= WHEN_TO_STOP) break;
  }
  if(y == 2*L+l)
  {
    if(VERBOSE) printf("OOPS, extended right all the way to %d\n",y);
    if(VERBOSE) fprintf(stderr,"OOPS, extended right all the way to %d\n",y);
    /* exit(1); */
  }
  y = besty;
  totalbestscore = besttotalbestscore;
  nrepeatocc = bestnrepeatocc;
  nactiverepeatocc = bestnactiverepeatocc;
  if(VERBOSE >= 2) printf("AFTER EXTENDING RIGHT, FAMILY %d IS:\n",R);
  if(VERBOSE >= 2) print_totalbestscore_right(y);
}

void extend_left()
{
  int w, n, bestscore, tempscore, offset;
  char a, besta;
  int newtotalbestscore, newtotalbestscore_a;
  besta = 0;

  /* initialize score[L%2][][] */
  w = L;
  for(n=0; n<N; n++)
  {
    bestscore = savebestscore[n] + CAPPENALTY;
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      if(score_of_besty[n][offset+MAXOFFSET] > bestscore)
        bestscore = score_of_besty[n][offset+MAXOFFSET];
    }
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      score[w%2][n][offset+MAXOFFSET] = bestscore; 
      if(offset < 0) score[w%2][n][offset+MAXOFFSET] += -offset * GAP;
      if(offset > 0) score[w%2][n][offset+MAXOFFSET] += offset * GAP;
    }
    bestbestscore[n] = bestscore;
  }

  /* extend left */
  bestw = L;
  compute_totalbestscore_left(w); 
  // RMH:
  if(VERBOSE >= 4) print_totalbestscore_left(w); 
  for(w=L-1; w>=0; w--) /* determine master[w] and score[w%2][][] */
  {
    newtotalbestscore = 0;
    for(a=0; a<4; a++) 
    {
      newtotalbestscore_a = 0;
      for(n=0; n<N; n++)
      {
        bestscore = bestbestscore[n] + CAPPENALTY;
        if(bestscore < 0) bestscore = 0;
        for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
        {
          tempscore = compute_score_left(w,n,offset,a);
          if(tempscore > bestscore)
            bestscore = tempscore;
        }
        newtotalbestscore_a += bestscore;
      }
      if(newtotalbestscore_a > newtotalbestscore)
      {
        newtotalbestscore = newtotalbestscore_a;
        besta = a;
      }
    }

    master[w] = besta;
    for(n=0; n<N; n++)
    {
      for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
        score[w%2][n][offset+MAXOFFSET] = compute_score_left(w,n,offset,besta);
    }

    compute_totalbestscore_left(w); 
    if(VERBOSE >= 4) print_totalbestscore_left(w); 
  
    if(w - bestw <= -WHEN_TO_STOP) break;
  }
  if(w == -1)
  {
    if(VERBOSE) printf("OOPS, extended left all the way to %d\n",w);
    if(VERBOSE) fprintf(stderr,"OOPS, extended left all the way to %d\n",w);
    /* exit(1); */
  }
  w = bestw;
  totalbestscore = besttotalbestscore;
  nrepeatocc = bestnrepeatocc;
  nactiverepeatocc = bestnactiverepeatocc;
  if(VERBOSE >= 2) printf("AFTER EXTENDING LEFT, FAMILY %d IS:\n",R);
  if(VERBOSE >= 2) print_totalbestscore_left(w);
  if(VERBOSE >= 2) fprintf(stderr,"AFTER EXTENDING LEFT, FAMILY %d IS:\n",R);
  if(VERBOSE >= 2) fprint_totalbestscore_left(w);
  
  masterstart[R] = bestw;
  masterend[R] = besty;
  for(w=bestw; w<=besty; w++) masters[R][w] = master[w];
}

void mask_headptr(struct llist **headptr) /* removed[x]=1 means l-mer x was removed */
{
  int x, y, smallh, h, h2, flip;
  int startmask, endmask;
  struct repeatllist **repeatheadptr;
  struct repeatllist *repeattmp, *nextrepeattmp;
  struct llist *tmp, *tmp2;
  struct posllist *postmp;
  int rrev = 0;

  /* build repeatheadptr */
  if((repeatheadptr = (struct repeatllist **) 
        malloc(SMALLHASH_SIZE*sizeof(*repeatheadptr))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  }
  // initialize repeatheadptr 
  for(smallh=0; smallh<SMALLHASH_SIZE; smallh++)
    repeatheadptr[smallh] = NULL;
  
  // Build a ds with all lmer's from one strand of the new master consensus R
  for(x=masterstart[R-1]; x<=masterend[R-1]-l+1; x++) /* l-mer masters[R-1]+x */
  {
    smallh = smallhash_function(masters[R-1]+x);
    if(smallh<0) continue;
    repeattmp = repeatheadptr[smallh];
    while(repeattmp != NULL)
    {
      /* already in there */
      if(lmermatch(repeattmp->value,masters[R-1]+x) == 1) break; 
      repeattmp = repeattmp->next;
    }
    /* add new guys */
    if(repeattmp != NULL) continue; /* already in there */
    if((repeattmp = (struct repeatllist *) malloc(sizeof(*repeattmp))) == NULL)
    {
      fprintf(stderr,"Out of memory\n");  exit(1);
    }
    for(y=0; y<l; y++) repeattmp->value[y] = masters[R-1][x+y];
    repeattmp->repeatpos = x;
    repeattmp->next = repeatheadptr[smallh]; /* either NULL, or l-mer with hash smallh */
    repeatheadptr[smallh] = repeattmp; 
  }

  /* find and extend hits */
  for(smallh=0; smallh<SMALLHASH_SIZE; smallh++)
  {
    repeattmp = repeatheadptr[smallh];
    while(repeattmp != NULL)
    {
      /* find tmp (in headptr) corresponding to repeattmp (in repeatheadptr) */
      h = hash_function(repeattmp->value);
      tmp = headptr[h];
      while(tmp != NULL)
      {
        // RMH: TODO: Possibly no longer needed
        if(lmermatch(sequence+(tmp->lastocc),repeattmp->value))
	{
          flip = 0;
          break;
        }
        else if(lmermatchrc(sequence+(tmp->lastocc),repeattmp->value))
	{
          flip = 1;
          break;
	}
        tmp = tmp->next;
      }
      if(tmp == NULL) { repeattmp = repeattmp->next; continue; } /* not in headptr */

      /* loop through pos of tmp: if not already masked, then extend and mask */ 
      postmp = tmp->pos;
      while(postmp != NULL)
      {
        // RMH: No longer needed because we are only keeping one strand
        //if(flip == 0) x = postmp->this;
        //else x = length-l - postmp->this; /* NOT length-1 */

        x = postmp->this;

        /* skip if already masked */
        if(removed[x] == 1) { postmp = postmp->next; continue; }
    
        // RMH: Fwd/Rev Mods
        rrev = 1;
        if(lmermatch(sequence+x,repeattmp->value))
          rrev = 0;

        /* extend this hit */
        // RMH: TODO: Fix lmer 2 or more times in repeat problem!!!
        //      This problem occurs when an lmer occurs 2 or more times
        //      in the model.  In these cases the model and the sequence
        //      lmer may not be matched up correctly.
        /* bestrepeaty in [0,2L], bestsequencey in [0,length] */
        bestsequencey = maskextend_right( rrev, x, repeattmp->repeatpos ); 
        /* bestrepeatw in [0,2L], bestsequencew in [0,length] */
        bestsequencew = maskextend_left( rrev, x, repeattmp->repeatpos );  

        /* mask the extended hit AND its rc */
        if(EXTRAMASK)
        {
          startmask = bestsequencew-l+1; 
          if(startmask<0) startmask = 0;
          endmask = bestsequencey;
        }
        else
        {
          startmask = bestsequencew;
          endmask = bestsequencey-l+1;
        }
        for(y=startmask; y<=endmask; y++) 
        {
          if(removed[y] == 1) continue; /* already removed */
          h2 = hash_function(sequence+y);
          if(h2<0) continue;
          tmp2 = headptr[h2];
          while(tmp2 != NULL)
          {
            /* forward or rc match */
            if(lmermatcheither(sequence+(tmp2->lastocc),sequence+y) == 1) 
            {      
              tmp2->freq -= 1;
              break;
            }
            tmp2 = tmp2->next;
          }
          removed[y] = 1;
          // RMH: No longer needed, Fwd/Rev Mods
          //removed[length-l-y] = 1;          
	}
        if(VERBOSE >= 3)
        {
          printf("Removed %d to %d (and rc):\n",startmask,endmask);
          for(y=startmask; y<=endmask; y++) 
              printf("%c",num_to_char(sequence[y]));
          printf("\n");
        }

        postmp = postmp->next;
      }
      repeattmp = repeattmp->next;
    }
  }
   
  /* free repeatheadptr */
  for(smallh=0; smallh<SMALLHASH_SIZE; smallh++)
  {
    repeattmp = repeatheadptr[smallh];
    while(repeattmp != NULL)
    {
      nextrepeattmp = repeattmp->next;
      free(repeattmp);
      repeattmp = nextrepeattmp;
    }
  }
  free(repeatheadptr);
}

// mask extend sequence right.  The model may be flipped either way 
// depending on the orientation of the sequence lmer.
// 
// fwd:   masterstart[R-1] modelpos      masterend[R-1]
//                |          |               |
//                v          v               v
//     Model    .............ACCGTAA...............
//     Seq   ................ACCGTAA.......................
//                           ^     *         ^
//                           |               |
//                         seqpos    seqpos + (masterend[R-1] - modelpos)
//
//
// rev:   masterstart[R-1] modelpos      masterend[R-1]
//                |          |               |
//                v          v               v
//     Model    .............ACCGTAA...............
//     Seq   ................TTACGGT.......................
//                           ^     *         ^
//                           |               |
//                       seqpos    seqpos + l + (modelpos - masterstart[R-1])
int maskextend_right( char rc, int seqpos, int modelpos )
{
  int offset, bestmaskscore, bestbestmaskscore, bestoffset, bestbestoffset;
  int bEnd, x, bestx;
  int maxext;
  char mbase;

  //printf("maskextend_right( %d, %d, %d ): Called\n", rc, seqpos, modelpos );

  for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
  {
    maskscore[1][offset+MAXOFFSET] = 0;  /* actually l*MATCH -- but arbitrary */
    if(offset < 0) maskscore[1][offset+MAXOFFSET] += -offset * GAP;
    if(offset > 0) maskscore[1][offset+MAXOFFSET] += offset * GAP;
  }

  bestbestmaskscore = 0; /* actually l*MATCH -- but arbitrary */ 
  bestbestoffset = 0;
  bestoffset = 0;
  bestx = 0; 

 
  // RMH: Calculate boundary start/end positions ( start inclusive, 
  //      end exclusive )
  bEnd = 0;
  x = 0;
  while ( boundaries[x] > 0  )
  {
   if ( seqpos < boundaries[x] + PADLENGTH )
   {
     bEnd = boundaries[x] + PADLENGTH;
     break;
   }
   x++;
  }

  maxext = masterend[R-1] - (modelpos + l + 1);
  if ( rc )
    maxext = modelpos - (masterstart[R-1] + 1);

  //printf("maskextend_right: bEnd = %d, maxext = %d, bestx = %d\n", bEnd, maxext, bestx );
  for( x = 0; x < maxext; x++ )
  {
    // RMH: Short circuit search at the end of the sequence boundary
    if ( seqpos + l + x == bEnd )
      break;

    bestmaskscore = -1000000000;

    if ( rc ) 
      mbase = compl(masters[R-1][modelpos - x - 1]);
    else
      mbase = masters[R-1][modelpos + l + x];

    //printf("mbase = %c/%c bestbestmaskscore = %d\n", num_to_char( mbase ), num_to_char( sequence[(seqpos+l+x)] ), bestbestmaskscore );
 
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      // RMH: Use sequence boundaries to limit extension
      if ( seqpos + l + x + offset >= bEnd )
      {
        if ( VERBOSE > 2 ) printf("Short circuiting mask ( %d ) due to "
                                  "sequence boundary ( bEnd = %d )\n", 
                         (seqpos+l+x+offset), bEnd );
        maskscore[(x+2)%2][offset+MAXOFFSET] = 0; 
      }else {
        maskscore[(x+2)%2][offset+MAXOFFSET] = 
                        compute_maskscore_right(mbase,x,seqpos,offset); 
      }

      if(maskscore[(x+2)%2][offset+MAXOFFSET] > bestmaskscore)
      {
        bestmaskscore = maskscore[(x+2)%2][offset+MAXOFFSET];
        bestoffset = offset;
      }
    }
    if(bestmaskscore > bestbestmaskscore)
    {
      bestbestmaskscore = bestmaskscore;
      bestbestoffset = bestoffset;
      bestx = x;
    }
    if(x - bestx >= WHEN_TO_STOP) break;
  }
  return ( seqpos + l + bestx + bestbestoffset );
}


// 
// fwd:   masterstart[R-1] modelpos      masterend[R-1]
//                |          |               |
//                v          v               v
//     Model    .............ACCGTAA...............
//     Seq   ................ACCGTAA.......................
//                ^          ^     
//                |          |      
//                |        seqpos   
//           seqpos - ( modelpos - masterstart[R-1] )
//
//
// rev:   masterstart[R-1] modelpos      masterend[R-1]
//                |          |               |
//                v          v               v
//     Model    .............ACCGTAA...............
//     Seq   ................TTACGGT.......................
//                ^          ^ 
//                |          | 
//                |        seqpos 
//             seqpos - ( masterend[R-1] - modelpos + l - 1)
int maskextend_left( char rc, int seqpos, int modelpos ) 
{
  int offset, bestmaskscore, bestbestmaskscore, bestoffset, bestbestoffset;
  int maxext, bStart, x, bestx;
  char mbase;

  //printf("maskextend_left( %d, %d, %d ): Called\n", rc, seqpos, modelpos );

  for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
  {
    maskscore[1][offset+MAXOFFSET] = 0;  /* actually saved best score -- but arbitrary */
    if(offset < 0) maskscore[1][offset+MAXOFFSET] += -offset * GAP;
    if(offset > 0) maskscore[1][offset+MAXOFFSET] += offset * GAP;
  }

  bestbestoffset = 0;
  bestoffset = 0;  /* Remove a compiler warning. */
  bestbestmaskscore = 0; /* actually saved best score -- but arbitrary */ 
  bestx = 0;

  bStart = 0;
  x = 0;
  while ( boundaries[x] > 0  )
  {
   if ( seqpos < boundaries[x] + PADLENGTH )
   {
     if ( x > 0 )
       bStart = boundaries[x-1] + PADLENGTH;
     break;
   }
   x++;
  }

  maxext = modelpos - masterstart[R-1];
  if ( rc ) // extending right down model to masterend
    maxext = masterend[R-1] - (modelpos + l - 1);
   
  //printf("maskextend_left: bStart = %d, maxext = %d\n", bStart, maxext );
  for(x = 0; x < maxext; x++ )
  {
    // RMH: Short circuit search at the edge of the sequence
    if ( seqpos - x - 1 == bStart )
       break;

    bestmaskscore = -1000000000;

    if ( rc ) 
      mbase = compl(masters[R-1][modelpos+l+x]);
    else
      mbase = masters[R-1][modelpos-x-1];

    //printf("mbase = %c/%c bestbestmaskscore = %d x = %d maxext = %d masterend[] = %d\n", num_to_char( mbase ), num_to_char( sequence[(seqpos-x-1)] ), bestbestmaskscore, x, maxext,  masterend[R-1] );

    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {

      // RMH: Use sequence boundaries to limit extension
      if ( seqpos-x-1+offset < bStart )
      {
        if ( VERBOSE > 2 ) printf("Short circuiting mask ( %d ) due to "
                                  "sequence boundary ( bStart=%d )\n", 
                          (seqpos-x-1+offset), bStart  );
        maskscore[(x+2)%2][offset+MAXOFFSET] = 0;
      }else {
        maskscore[(x+2)%2][offset+MAXOFFSET] = compute_maskscore_left(mbase,x,seqpos,offset); 
      }

      if(maskscore[(x+2)%2][offset+MAXOFFSET] > bestmaskscore)
      {
        bestmaskscore = maskscore[(x+2)%2][offset+MAXOFFSET];
        bestoffset = offset;
      }
    }
    if(bestmaskscore > bestbestmaskscore)
    {
      bestbestmaskscore = bestmaskscore;
      bestbestoffset = bestoffset;
      bestx = x;
    }
    if( x - bestx >= WHEN_TO_STOP ) break; 
  }
  return( seqpos - bestx + bestbestoffset - 1);
}

/* ************************************************************************ */

int hash_function2(char *lmer) /* IS symmetric wrt reverse complements */
{
  int x, ans, ans2;

  for(x=0; x<l; x++)
  {
    if(char_to_num(lmer[x]) == 99) return -1;
  }
  ans = 0;
  for(x=0; x<l; x++)
    ans = (4*ans + (char_to_num(lmer[x])%4)) % HASH_SIZE;
  ans2 = 0;
  for(x=0; x<l; x++)
    ans2 = (4*ans2 + ((3-char_to_num(lmer[l-1-x]))%4)) % HASH_SIZE;
  if(ans2 > ans) ans = ans2;
  return ans;
}

int hash_function(char *lmer) /* IS symmetric wrt reverse complements */
{
  int x, ans, ans2;

  for(x=0; x<l; x++)
  {
    if(lmer[x] == 99) return -1;
  }
  ans = 0;
  for(x=0; x<l; x++)
    ans = (4*ans + (lmer[x]%4)) % HASH_SIZE;
  ans2 = 0;
  for(x=0; x<l; x++)
    ans2 = (4*ans2 + ((3-lmer[l-1-x])%4)) % HASH_SIZE;
  if(ans2 > ans) ans = ans2;
  return ans;
}

int smallhash_function(char *lmer) /* IS symmetric wrt reverse complements */
{
  int x, ans;

  for(x=0; x<smalll; x++)
  {
    if(lmer[x] == 99) return -1;
  }
  ans = lmer[0];
  for(x=1; x<smalll; x++)
    ans = (4*ans + (lmer[x]%4)) % SMALLHASH_SIZE;
  return ans;
}

int lmermatch(char *lmer1, char *lmer2)
{
  int x;

  for(x=0; x<l; x++)
  {
    if(lmer1[x] != lmer2[x]) return 0;
  }
  return 1;
}

int lmermatcheither(char *lmer1, char *lmer2) /* forward or rc match */
{
  int x;

  int lmermatchrc(char *lmer1, char *lmer2);

  for(x=0; x<l; x++)
  {
    if(lmer1[x] != lmer2[x]) return lmermatchrc(lmer1, lmer2);
  }
  return 1;
}

int lmermatchrc(char *lmer1, char *lmer2) /* rc match */
{
  int x;

  for(x=0; x<l; x++)
  {
    if((lmer1[x] + lmer2[l-1-x] != 3)) return 0;
  }
  return 1;
}

int mismatches(char *lmer1, char *lmer2)
{
  int x, ans;

  ans = 0;
  for(x=0; x<l; x++)
  {
    if(lmer1[x] != lmer2[x]) ans++;
  }
  return ans;
}

int is_degenerate(char *this)
{
  /* defn of degenerate: any repeated pattern of length 2-5 (see p. 157) */
  /* ALSO: need all 4 lett */

  int i; /* length of repeated pattern */
  int x, undeg;
  int occ4[4];
  int MISMATCHES;

  MISMATCHES = 5;

  /* first, check all 4 lett */
  for(x=0; x<4; x++) occ4[x] = 0;
  for(x=0; x<l; x++) occ4[ (int)this[x] ] = 1;
  for(x=0; x<4; x++)
  {
    if(occ4[x] == 0) return 1;
  }

  for(i=2; i<=5; i++) /* from beg */
  {
    undeg = 0;
    for(x=i; x<l; x++)
    {
      if(this[x] != this[(x%i)])
        undeg++;
    }
    if(undeg <= MISMATCHES) return 1;
  }
  for(i=2; i<=5; i++) /* from end */
  {
    undeg = 0;
    for(x=0; x<l-i; x++)
    {
      if(this[x] != this[(x%i)+l-i])
        undeg++;
    }
    if(undeg <= MISMATCHES) return 1;
  }
  return 0;
}

void compute_totalbestscore_right(int y)
{
  int n, bestscore, offset;

  nrepeatocc = 0;
  nactiverepeatocc = 0;
  totalbestscore = 0;
  for(n=0; n<N; n++)
  {
    bestscore = bestbestscore[n] + CAPPENALTY;
    if(bestscore < 0) bestscore = 0;
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      if(score[y%2][n][offset+MAXOFFSET] > bestscore)
        bestscore = score[y%2][n][offset+MAXOFFSET];
    }
    if(bestscore > 0) nrepeatocc++;
    if(bestscore > bestbestscore[n] + CAPPENALTY) nactiverepeatocc++;
    if(bestscore > bestbestscore[n]) bestbestscore[n] = bestscore;
    totalbestscore += bestscore;
  }

  if((totalbestscore >= besttotalbestscore + (y-besty)*MINIMPROVEMENT)
     && (totalbestscore > besttotalbestscore)) /* to deal with MINIMPROVEMENT=0 */
  {
    besty = y;
    besttotalbestscore = totalbestscore;
    bestnrepeatocc = nrepeatocc;
    bestnactiverepeatocc = nactiverepeatocc;
    for(n=0; n<N; n++) 
    {
      savebestscore[n] = bestbestscore[n];
      for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
        score_of_besty[n][offset+MAXOFFSET] = score[y%2][n][offset+MAXOFFSET];
    }
  }

  return;
}

void print_totalbestscore_right(int y)
{
  int x;
  char num_to_char(char z);

  printf("%d letters total:\n",y+1-L);
  for(x=L; x<=y; x++) 
  {
    printf("%c",num_to_char(master[x]));
    if((x-L)%80 == 79) printf("\n");
  }
  printf("\n");
  printf("totalbestscore is %d (%d of %d are good, %d of %d are active)\n\n",totalbestscore,nrepeatocc,N,nactiverepeatocc,N);

  return;
}

void compute_totalbestscore_left(int w)
{
  int n, bestscore, offset;

  nrepeatocc = 0;
  nactiverepeatocc = 0;
  totalbestscore = 0;
  for(n=0; n<N; n++)
  {
    bestscore = bestbestscore[n] + CAPPENALTY;
    if(bestscore < 0) bestscore = 0;
    for(offset=-MAXOFFSET; offset<=MAXOFFSET; offset++)
    {
      if(score[w%2][n][offset+MAXOFFSET] > bestscore)
        bestscore = score[w%2][n][offset+MAXOFFSET];
    }
    if(bestscore > 0) nrepeatocc++;
    if(bestscore > bestbestscore[n] + CAPPENALTY) nactiverepeatocc++;
    if(bestscore > bestbestscore[n]) bestbestscore[n] = bestscore;
    totalbestscore += bestscore;
  }

  if((totalbestscore >= besttotalbestscore + (bestw-w)*MINIMPROVEMENT)
     && (totalbestscore > besttotalbestscore)) /* to deal with MINIMPROVEMENT=0 */
  {
    bestw = w;
    besttotalbestscore = totalbestscore;
    bestnrepeatocc = nrepeatocc;
    bestnactiverepeatocc = nactiverepeatocc;
  }

  return;
}

void print_totalbestscore_left(int w)
{
  int x;
  char num_to_char(char z);

  printf("%d letters total:\n",besty-w+1);
  for(x=w; x<=besty; x++) 
  {
    printf("%c",num_to_char(master[x]));
    if((x-w)%80 == 79) printf("\n");
  }
  printf("\n");
  printf("totalbestscore is %d (%d of %d are good, %d of %d are active)\n\n",totalbestscore,nrepeatocc,N,nactiverepeatocc,N);

  return;
}

void fprint_totalbestscore_left(int w)
{
  int x;
  char num_to_char(char z);

  fprintf(stderr,"%d letters total:\n",besty-w+1);
  for(x=w; x<=besty; x++) 
  {
    fprintf(stderr,"%c",num_to_char(master[x]));
    if((x-w)%80 == 79) fprintf(stderr,"\n");
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"totalbestscore is %d (%d of %d are good, %d of %d are active)\n\n",totalbestscore,nrepeatocc,N,nactiverepeatocc,N);

  return;
}

int compute_score_right(int y, int n, int offset, char a)
{
  int oldoffset, tempscore, ans, ismatch, x;
  int bStart, bEnd;

  //RMH: Sequence boundaries check.  
  //      bStart                  bEnd
  //    inclusive---------------exclusive
  bStart = PADLENGTH;
  if ( upperBoundI[n] > 0 )
    bStart = boundaries[upperBoundI[n]-1] + PADLENGTH;
  bEnd = boundaries[upperBoundI[n]] + PADLENGTH;
  bStart =0;
  bEnd = 50000000;

  if ( rev[n] ) // going towards bStart
  {
    if( pos[n]-(offset+y-L-l)-1 < bStart  )
    {
      if ( VERBOSE > 2 ) printf("Short circuiting lmer %d  ( %d ) due to sequence boundary ( %d<- ) upperBoundI = %d \n", n, (pos[n]-(offset+y-L-l)-1), bStart, upperBoundI[n] );
      return 0;
    }
  }
  else // going towards bEnd
  {
    if( pos[n]+offset+y-L >= bEnd )
    {
      if ( VERBOSE > 2 ) printf("Short circuiting lmer %d due to sequence boundary ( ->%d )\n", n, bEnd );
      return 0;
    }
  }

  ans = -1000000000; 
  /* first, try oldoffset = offset+1.  
     This means master[y]=a aligns to - */
  if(offset < MAXOFFSET)
  {
    oldoffset = offset+1;
    tempscore = score[(y-1)%2][n][oldoffset+MAXOFFSET] + GAP;
    if(tempscore > ans) ans = tempscore;
  }

  /* next, try oldoffset = offset.  
     This means master[y]=a aligns to sequence[pos[n]+offset+y-L] */
  oldoffset = offset;
  tempscore = score[(y-1)%2][n][oldoffset+MAXOFFSET];
  if ( rev[n] )
  {
    if(a == compl(sequence[pos[n]-(offset+y-L-l)-1]))
    {
      tempscore += MATCH;
    }
    else
    {
      tempscore += MISMATCH;
    }
  }
  else 
  {
    if(a == sequence[pos[n]+offset+y-L])
      tempscore += MATCH;
    else
      tempscore += MISMATCH;
  }
  if(tempscore > ans) ans = tempscore;

  /* finally, try oldoffset < offset 
    This means master[y]=a aligns to one of 
    sequence[pos[n]+oldoffset+y-L] ... sequence[pos[n]+offset+y-L]. 
    Also, offset-oldoffset gaps */
  for(oldoffset=-MAXOFFSET; oldoffset<offset; oldoffset++)
  {
    ismatch = 0;
    for(x=oldoffset; x<=offset; x++)
    {
      if ( rev[n] )
      {
        if(a == compl(sequence[pos[n]-(x+y-L-l)-1])) ismatch = 1;
      }
      else
      {
        if(a == sequence[pos[n]+x+y-L]) ismatch = 1;
      }
    }
    tempscore = score[(y-1)%2][n][oldoffset+MAXOFFSET];
    tempscore += (offset - oldoffset) * GAP;
    if(ismatch == 1) tempscore += MATCH;
    else tempscore += MISMATCH;
    if(tempscore > ans) ans = tempscore;
  }

  return ans;
}

//
// RMH: Modified to handle reverse strand lmers.  The rev[n] structure
//      is global and holds a flag indicating if pos[n] points to a 
//      forward or reverse strand lmer.
//
int compute_score_left(int w, int n, int offset, char a )
{
  int oldoffset, tempscore, ans, ismatch, x;
  int bStart, bEnd;

  //RMH: Sequence boundaries check.  If bEnd is set
  //      then use boundaries to short circuit calc and return 0.
  //      bStart                  bEnd
  //    inclusive---------------exclusive
  bStart = PADLENGTH;
  if ( upperBoundI[n] > 0 )
    bStart = boundaries[upperBoundI[n]-1] + PADLENGTH;
  bEnd = boundaries[upperBoundI[n]] + PADLENGTH;
  bStart =0;
  bEnd = 50000000;
  if ( rev[n] ) // going towards bEnd
  {
    if( pos[n]-(offset+w-L-l)-1 >= bEnd  )
    {
      if ( VERBOSE > 2 ) printf("Short circuiting lmer %d  ( %d ) due to "
                                "sequence boundary ( ->%d ) upperBoundI "
                                "= %d \n", n, (pos[n]-(offset+w-L-l)-1), 
                                bEnd, upperBoundI[n] );
      return 0;
    }
  }
  else // going towards bStart
  {
    if( pos[n]+offset+w-L < bStart )
    {
      if ( VERBOSE > 2 ) printf("Short circuiting lmer %d due to sequence"
                                " boundary ( %d<- )\n", n, bStart );
      return 0;
    }
  }

  ans = -1000000000; 
  /* first, try oldoffset = offset-1.  
     This means master[w]=a aligns to - */
  if(offset > -MAXOFFSET)
  {
    oldoffset = offset-1;
    tempscore = score[(w+1)%2][n][oldoffset+MAXOFFSET] + GAP;
    if(tempscore > ans) ans = tempscore;
  }

  /* next, try oldoffset = offset.  
     This means master[w]=a aligns to sequence[pos[n]+offset+w-L] */
  oldoffset = offset;
  tempscore = score[(w+1)%2][n][oldoffset+MAXOFFSET];

  if ( rev[n] )
  {
    if(a == compl(sequence[pos[n]-(offset+w-L-l)-1]))
      tempscore += MATCH;
    else
      tempscore += MISMATCH;
  }
  else 
  {
    if(a == sequence[pos[n]+offset+w-L])
      tempscore += MATCH;
    else
      tempscore += MISMATCH;
  }
  if(tempscore > ans) ans = tempscore;

  /* finally, try oldoffset > offset 
    This means master[w]=a aligns to one of 
    sequence[pos[n]+offset+w-L] ... sequence[pos[n]+oldoffset+w-L]. 
    Also, oldoffset-offset gaps */
  for(oldoffset=offset+1; oldoffset<=MAXOFFSET; oldoffset++)
  {
    ismatch = 0;
    for(x=offset; x<=oldoffset; x++)
    {
      if ( rev[n] )
      {
        if(a == compl(sequence[pos[n]-(x+w-L-l)-1])) ismatch = 1;
      }
      else
      {
        if(a == sequence[pos[n]+x+w-L]) ismatch = 1;
      }
    }
    tempscore = score[(w+1)%2][n][oldoffset+MAXOFFSET];
    tempscore += (oldoffset - offset) * GAP;
    if(ismatch == 1) tempscore += MATCH;
    else tempscore += MISMATCH;
    if(tempscore > ans) ans = tempscore;
  }

  return ans;
}

int compute_maskscore_right(char mbase, int repeaty, int seqpos, int offset)
{
  int oldoffset, tempscore, ans, ismatch, x;

  ans = -1000000000; 
  /* first, try oldoffset = offset+1.  
     This means masters[R-1][repeaty] (mbase) aligns to - */
  if(offset < MAXOFFSET)
  {
    oldoffset = offset+1;
    tempscore = maskscore[(repeaty+1)%2][oldoffset+MAXOFFSET] + GAP;
    if(tempscore > ans) ans = tempscore;
  }

  /* next, try oldoffset = offset.  
     This means masters[R-1][repeaty] (mbase) aligns to sequence[repeaty+repeat2sequence+offset] */
  oldoffset = offset;
  tempscore = maskscore[(repeaty+1)%2][oldoffset+MAXOFFSET];
  if(mbase == sequence[seqpos + l + repeaty + offset])
    tempscore += MATCH;
  else
    tempscore += MISMATCH;
  if(tempscore > ans) ans = tempscore;

  /* finally, try oldoffset < offset 
    This means masters[R-1][repeaty] (mbase) aligns to one of 
    sequence[repeaty+repeat2sequence+oldoffset] ... sequence[repeaty+repeat2sequence+offset]. 
    Also, offset-oldoffset gaps */
  for(oldoffset=-MAXOFFSET; oldoffset<offset; oldoffset++)
  {
    ismatch = 0;
    for(x=oldoffset; x<=offset; x++)
    {
      if(mbase == sequence[seqpos + l + repeaty + x]) ismatch = 1;
    }
    tempscore = maskscore[(repeaty+1)%2][oldoffset+MAXOFFSET];
    tempscore += (offset - oldoffset) * GAP;
    if(ismatch == 1) tempscore += MATCH;
    else tempscore += MISMATCH;
    if(tempscore > ans) ans = tempscore;
  }

  return ans;
}

int compute_maskscore_left(char mbase, int repeatw, int seqpos, int offset)
{
  int oldoffset, tempscore, ans, ismatch, x;

  ans = -1000000000; 
  /* first, try oldoffset = offset-1.  
     This means masters[R-1][repeatw] aligns to - */
  if(offset > -MAXOFFSET)
  {
    oldoffset = offset-1;
    tempscore = maskscore[(repeatw+1)%2][oldoffset+MAXOFFSET] + GAP;
    if(tempscore > ans) ans = tempscore;
  }
  //printf("off = %d, ans = %d, ", offset, ans );

  /* next, try oldoffset = offset.  
     This means masters[R-1][repeatw] aligns to sequence[repeatw+repeat2sequence+offset] */
  oldoffset = offset;
  tempscore = maskscore[(repeatw+1)%2][oldoffset+MAXOFFSET];
  if(mbase == sequence[seqpos-repeatw-1+offset])
    tempscore += MATCH;
  else
    tempscore += MISMATCH;
  if(tempscore > ans) ans = tempscore;
  //printf("ans = %d\n", ans );

  /* finally, try oldoffset > offset 
    This means masters[R-1][repeatw] aligns to one of 
    sequence[repeatw+repeat2sequence+offset] ... sequence[repeatw+repeat2sequence+oldoffset]. 
    Also, oldoffset-offset gaps */
  for(oldoffset=offset+1; oldoffset<=MAXOFFSET; oldoffset++)
  {
    ismatch = 0;
    for(x=offset; x<=oldoffset; x++)
    {
      if(mbase == sequence[seqpos-repeatw-1+x]) ismatch = 1;
    }
    tempscore = maskscore[(repeatw+1)%2][oldoffset+MAXOFFSET];
    tempscore += (oldoffset - offset) * GAP;
    if(ismatch == 1) tempscore += MATCH;
    else tempscore += MISMATCH;
    if(tempscore > ans) ans = tempscore;
  }

  return ans;
}

int minimum_edit_distance(char *lmer, int r)
{
  int x, dist, mindist;
  char lmerrc[l];

  int edit_distance(char *lmer1, char *lmer2);

  for(x=0; x<l; x++) lmerrc[x] = 3 - lmer[l-1-x];
  mindist = 1000;
  for(x=masterstart[r]; x<=masterend[r]-l+1; x++)
  {
    dist = edit_distance(lmer, masters[r]+x);
    if(dist < mindist) mindist = dist;
    dist = edit_distance(lmerrc, masters[r]+x);
    if(dist < mindist) mindist = dist;
  }

  return mindist;
}

int edit_distance(char *lmer1, char *lmer2)
{
  int x, y, thisdist, PARTIAL_EDIT_OK;

  x=0; for(y=0; y<=l; y++) distmat[x][y] = y;
  for(x=1; x<=l; x++) /* row x */
  {
    distmat[x][0] = x;
    for(y=1; y<=l; y++)
    {
      thisdist = distmat[x-1][y-1]; if(lmer1[x-1] != lmer2[y-1]) thisdist++;
      if(distmat[x-1][y]+1 < thisdist) thisdist = distmat[x-1][y]+1;
      if(distmat[x][y-1]+1 < thisdist) thisdist = distmat[x][y-1]+1;
      distmat[x][y] = thisdist;
    }
  }
  thisdist = distmat[l][l];
  PARTIAL_EDIT_OK = 1;
  if(PARTIAL_EDIT_OK) /* any entry on last row or column */
  {
    for(x=0; x<l; x++) { if(distmat[x][l] < thisdist) thisdist = distmat[x][l]; }
    for(y=0; y<l; y++) { if(distmat[l][y] < thisdist) thisdist = distmat[l][y]; }
  }
  return thisdist;
}

double compute_entropy(char *lmer)
{
  int x;
  int count[4];  
  double answer, y;
  
  for(x=0; x<4; x++) count[x] = 0;
  for(x=0; x<l; x++) count[ (int)lmer[x]] += 1;
  answer = 0.0;
  for(x=0; x<4; x++)
  {
    if(count[x] == 0) continue;
    y = ((double)count[x])/((double)l);
    answer += y*log(y);
  }
  return answer;
}

char char_to_num(char c)
{
  if(c == 'A') return 0;
  if(c == 'C') return 1;
  if(c == 'G') return 2;
  if(c == 'T') return 3;
  if(c == 'a') return 0;
  if(c == 'c') return 1;
  if(c == 'g') return 2;
  if(c == 't') return 3;
  if(c == 'N') return 99;
  if(IUPAC(c)) return 99;
  if(c == 'n') return 99;
  if(c == 'x') return 99;
  fprintf(stderr,"ERROR3: c=%d->%c\n",c,c);
  exit(1);
}

char num_to_char(char z)
{
  if(z == 0) return (char) 'A';
  if(z == 1) return (char) 'C';
  if(z == 2) return (char) 'G';
  if(z == 3) return (char) 'T';
  return (char) 'N';
}

// RMH: Fwd/Rev Mods
char compl( char c )
{
  if(c == 0) return 3;
  if(c == 1) return 2;
  if(c == 2) return 1;
  if(c == 3) return 0;
  return 99;
}

/* 
 * Set l = ceil( 1 + log(\sum(g.segments[*].length)/log(4) );
 */
int default_l(int len) {
	return ceil( 1 + log(len) / log(4) );
}
