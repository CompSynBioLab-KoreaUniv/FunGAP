#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "build_lmer_table.h"
#include "cmd_line_opts.h"

int default_l(int);

/* build_lmer_table.c -- build table of l-mers */

// RMH: From version.c which in turn is created by the Makefile
extern const char *Version;

#define PADLENGTH 11000 /* important -- must match build_repeat_families.c */
#define HASH_SIZE 16000057 /* prime */
#define BLOCKLENGTH 1000000000
#define IUPAC(c) c == 'R' || c=='r' || c=='Y' || c=='y' || c=='M' || c=='m' || c=='K' || c=='k' || c=='W' || c=='w' || c=='S' || c=='s' || c=='B' || c=='b' || c=='D' || c=='d' || c=='H' || c=='h' || c=='V' || c=='v'

int length;                 /* length of genome sequence */
int l;                      /* length of l-mers to look at.  Odd or even is Ok. */
char* LMER_TABLE_FILE;      /* Where to put the results... */
char* SEQUENCE_FILE;        /* Where to get the sequence from... */
int VERBOSE    = 0;         /* How chatty should I be?  Really chatty?  Really really chatty? */
int TANDEMDIST = 500;       /* Distance between countable instances of an l-mer. */
int SORT       = 0;
int  MAXLENGTH;             /* Length of sequence to allow. */
int MINTHRESH  = 3;
int START = -1;
int STOP  = -1;

char *sequence;

void usage() 
{
	fprintf(stderr, "build_lmer_table Version %s\n\n"\
		"Usage:\n"\
		"   build_lmer_table -l <l> -sequence <seq> -freq <output> [opts]\n"\
		"       -tandem <d> --- tandem distance window (def: 500)\n"\
		"       -min <#> --- smallest number of required lmers (def: 3)\n"\
		"       -v --- output a small amount of debugging information.\n",
                Version
	);

	exit(1);
}

int main(int argc, char* argv[])
{
  time_t start;
  struct llist **headptr;
  int *sortedocc, *sortedfreq, ngoodlmers;
  long *sortedindex;
  FILE* fp;

  start = time(0);

  /* Gather command-line options... */
  if( 0 == 1*
      co_get_string(argc, argv, "-sequence", &SEQUENCE_FILE) *
      co_get_string(argc, argv, "-freq",     &LMER_TABLE_FILE)
  ) {
      usage();
  }

  co_get_int(argc, argv, "-tandem", &TANDEMDIST) || (TANDEMDIST = 500);
  co_get_int(argc, argv, "-min", &MINTHRESH) || (MINTHRESH = 3);

  co_get_bool(argc, argv, "-sort", &SORT);
  co_get_bool(argc, argv, "-v", &VERBOSE);

  co_get_int(argc, argv, "-start", &START);
  co_get_int(argc, argv, "-stop", &STOP);

  fp = fopen( SEQUENCE_FILE, "ro" );
  if( NULL == fp ) {
     fprintf(stderr, "Could not open sequence file %s\n", SEQUENCE_FILE);
     exit(1);
  }
  fseek(fp, 0, SEEK_END);
  MAXLENGTH = ftell(fp);  /* This is an approximation, but an overestimate, so we're ok! */
  fclose(fp);

  co_get_int(argc, argv, "-l", &l) || (l = default_l(MAXLENGTH));

  sequence = (char *)malloc( (MAXLENGTH + PADLENGTH) * sizeof(char) );
  if( NULL == sequence ) {
     fprintf(stderr, "Unable to allocate %d bytes for the sequence\n", MAXLENGTH + PADLENGTH );
     exit(1);
  }

  /* build sequence */
  length = build_sequence(sequence,SEQUENCE_FILE);

  /* build headptr */
  if((headptr = (struct llist **) malloc(HASH_SIZE*sizeof(*headptr))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  if(VERBOSE) fprintf(stderr,"  Done allocating headptr\n");
  build_headptr(headptr); if(VERBOSE) fprintf(stderr,"  Done building headptr\n");

  /* sort headptr */
  ngoodlmers = count_headptr(headptr);
  if(VERBOSE) fprintf(stderr,"  There are %d l-mers\n",ngoodlmers);
  if((sortedfreq = (int *) malloc(ngoodlmers*sizeof(*sortedfreq))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  if((sortedocc = (int *) malloc(ngoodlmers*sizeof(*sortedocc))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  if((sortedindex = (long *) malloc(ngoodlmers*sizeof(*sortedindex))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 
  sort_headptr(headptr, sortedfreq, sortedocc, sortedindex, ngoodlmers);
  if(VERBOSE) fprintf(stderr,"  Done sorting headptr\n");
  if(ngoodlmers == 0) { fprintf(stderr,"OOPS no good lmers\n"); exit(1); }

  print_lmers(headptr, sortedfreq, sortedocc, sortedindex, ngoodlmers);

  exit(0);
}

/* ************************************************************** */

int build_sequence(char *sequence, char *filename)
{
  int i, j;
  char c;
  FILE *fp;

  if( (fp = fopen(filename, "r")) == NULL)
  {
    fprintf(stderr,"Could not open input file %s\n", filename);  exit(1);
  }
  for(i=0; i<PADLENGTH; i++) sequence[i] = 99;
  i = 0;
  while( !feof(fp) && (i<MAXLENGTH) )
  { /* process one line of file */
    c = getc(fp);
    if(c == EOF) continue;
    if(c == '\n') continue;
    if(c == '>')
    {
      for(j=0; (getc(fp) != '\n') && !feof(fp); j++)
        ;
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
  fclose(fp);
  return (i+PADLENGTH); /* length of genome */
}

void build_headptr(struct llist **headptr)
{
  int h, i;
  int beg = START;
  int end = STOP;
  struct llist *tmp;

  if( START == -1 ) {
	beg = 0;
  }

  if( STOP == -1 ) {
        end = length - l + 1;
  }

  /* set each head pointer *headptr to NULL */
  for(h=0; h<HASH_SIZE; h++)
    headptr[h] = NULL;

  for(i=beg; i<end; i++)
  {
    if((i%BLOCKLENGTH == 0) && (i>0)) trim_headptr(headptr);
    if( !(i%100000) ) {
	fprintf(stderr, "%d      \r", i);
    }

    /* processing the l-mer sequence[i], ..., sequence[i+l-1] */
    h = hash_function(sequence+i);
    if(h<0) continue;
    tmp = headptr[h];
    while(tmp != NULL)
    {
      if(lmermatch(sequence+(tmp->lastplusocc),sequence+i) == 1) /* forward match */
      {      
        if(i - tmp->lastplusocc >= TANDEMDIST) 
        {
          tmp->freq += 1;
	}
        tmp->lastplusocc = i;
        break;
      }
      else if(lmermatchrc(sequence+(tmp->lastplusocc),sequence+i) == 1) /* rc match */
      {      
        if(i - tmp->lastminusocc >= TANDEMDIST) 
        {
          tmp->freq += 1;
	}
        tmp->lastminusocc = i;
        break;
      }
      tmp = tmp->next;
    }
    /* add new guys */
    if(tmp != NULL) continue; /* found l-mer, don't need to add in */
    if((tmp = (struct llist *) malloc(sizeof(*tmp))) == NULL)
    {
      fprintf(stderr,"Out of memory at i=%d\n",i);  exit(1);
    }
    tmp->lastplusocc = i;
    tmp->lastminusocc = -1000000;
    tmp->freq = 1;
    tmp->next = headptr[h]; /* either NULL, or some other l-mer with hash h */
    headptr[h] = tmp; 
  }
  trim_headptr(headptr);
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

int count_headptr(struct llist **headptr)
{
  int n, h;
  struct llist *tmp;

  n = 0;
  for(h=0; h<HASH_SIZE; h++)
  {
    tmp = headptr[h];
    while(tmp != NULL)
    {
      n++;
      tmp = tmp->next;
    }
  }
  return n;
}

void sort_headptr(struct llist **headptr, int *sortedfreq, int *sortedocc, long *sortedindex, int ngoodlmers)
{
  int n, h;
  int *order;
  struct llist *tmp;

  void q_sort(int *sortedfreq, int *sortedocc, long *sortedindex, int *order, int left, int right);
  long compute_index (char *this);
  char num_to_char(char z);

  if((order = (int *) malloc(ngoodlmers*sizeof(*order))) == NULL)
  {
    fprintf(stderr,"Out of memory\n");  exit(1);
  } 

  n = 0;
  for(h=0; h<HASH_SIZE; h++)
  {
    tmp = headptr[h];
    while(tmp != NULL)
    {
      sortedocc[n] = tmp->lastplusocc;
      sortedfreq[n] = tmp->freq;
      sortedindex[n] = compute_index(sequence+tmp->lastplusocc);
      n++;
      tmp = tmp->next;
    }
  }

  for(n=0; n<ngoodlmers; n++) order[n] = n;
  if(SORT) {
	q_sort(sortedfreq, sortedocc, sortedindex, order, 0, ngoodlmers - 1);
  }
   
}

void print_lmers(struct llist **headptr, int *sortedfreq, int *sortedocc, long *sortedindex, int ngoodlmers)
{
  int n, x;
  char lmer[l];
  FILE *fp;

  if( (fp = fopen(LMER_TABLE_FILE, "w")) == NULL)
  {
    fprintf(stderr,"Could not open file %s\n", LMER_TABLE_FILE);  exit(1);
  }


  for(n=0; n<ngoodlmers; n++)
  {
    for(x=0; x<l; x++) lmer[x] = sequence[sortedocc[n]+x];
    reverse_if_necessary(lmer);
    for(x=0; x<l; x++) fprintf(fp,"%c",num_to_char(lmer[x]));
    fprintf(fp,"\t%d\t%d\n",sortedfreq[n],sortedocc[n]);
  }
  fclose(fp);
}

void q_sort(int *sortedfreq, int *sortedocc, long *sortedindex, int *order, int left, int right)
{
  int pivot, l_hold, r_hold;
  int origfreq, origorder, origocc;
  long origindex;

  l_hold = left;
  r_hold = right;
  origfreq = sortedfreq[left]; 
  origorder = order[left]; 
  origocc = sortedocc[left]; 
  origindex = sortedindex[left];
  while (left < right)
  {
    while ((sortedfreq[right] <= origfreq) && (left < right))
      right--;
    if (left != right)
    {
      sortedfreq[left] = sortedfreq[right]; 
      order[left] = order[right];
      sortedocc[left] = sortedocc[right];
      sortedindex[left] = sortedindex[right];
      left++;
    }
    while ((sortedfreq[left] >= origfreq) && (left < right))
      left++;
    if (left != right)
    {
      sortedfreq[right] = sortedfreq[left]; 
      order[right] = order[left];
      sortedocc[right] = sortedocc[left];
      sortedindex[right] = sortedindex[left];
      right--;
    }
  }
  sortedfreq[left] = origfreq;      
  order[left] = origorder;
  sortedocc[left] = origocc;
  sortedindex[left] = origindex;

  pivot= left;
  left = l_hold;
  right = r_hold;
  if (left < pivot)
    q_sort(sortedfreq, sortedocc, sortedindex, order, left, pivot-1);
  if (right > pivot)
    q_sort(sortedfreq, sortedocc, sortedindex, order, pivot+1, right);
}

/* ************************************************************************ */

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

int lmermatch(char *lmer1, char *lmer2) /* forward match */
{
  int x;

  for(x=0; x<l; x++)
  {
    if(lmer1[x] != lmer2[x]) return 0;
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

long compute_index (char *this)
{
  int x;
  long index;
  char lmer[l];

  for(x=0; x<l; x++) lmer[x] = this[x];
  reverse_if_necessary(lmer);
  index = 0;
  for(x=0; x<l; x++)
    index = 4*index + lmer[x];

  return index;
}

void reverse_if_necessary(char *this)
{
  int x;
  char rcthis[l];

  for(x=0; x<l; x++) rcthis[x] = 3 - this[l-1-x];
  for(x=0; x<l; x++)
  {
    if(this[x] < rcthis[x]) return;
    if(this[x] > rcthis[x])
    {
      for(x=0; x<l; x++) this[x] = rcthis[x];
      return;
    }
  }
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
  if(c == 'M') return 99;
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

/* 
 * Set l = ceil( 1 + log(\sum(g.segments[*].length)/log(4) );
 */
int default_l(int len) {
	return ceil( 1 + log(len) / log(4) );
}
