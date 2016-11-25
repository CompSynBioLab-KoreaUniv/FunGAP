
#include "msps.h"


int GetOutputInterval(char *, int, int);



int main (int argc, char *argv[]) {
  FILE *msp_file, *seq_list; /* input files */
  FILE **output, *err, *msp_no; /* output files */
  int noof=1, i; /* noof == Number Of Output Files */
  char line[150], output_name[50];
  MSP_t cur;
  int img_ct=-1;

  /* checking command line */
  if (argc < 3) {
    printf("usage:\n  imagespread seq_name_list msp_list number_of_output_files\nlast parameter optional with default of 1\n");
    exit(1);
  }

  /* checking input */
  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  msp_file = fopen(argv[2], "r");
  if (!msp_file) {
    printf("Input MSP file %s not found.  Exit.\n", argv[2]);
    exit(2);
  }
  if (argc == 4) {noof = atoi(argv[3]);}
  if (noof < 1) noof = 1;

  /*files for output */
  err = fopen("images/errors", "w");
  if (!err) {
    printf("Can not open images/errors for output.  Exit.\n");
    exit(3);
  }
  msp_no = fopen("summary/ori_msp_no", "w");
  if (!msp_no) {
    printf("Can not open images/msp_no for output.  Exit.\n");
    exit(3);
  }
  output = (FILE **) malloc(noof*sizeof(FILE *));
  for (i=0; i<noof; i++) {
    sprintf(output_name, "images/spread%d", i+1);
    *(output+i) = fopen(output_name, "w");
    if (!*(output+i)) {
      printf("Can not open %s for output.  Exit.\n", output_name);
      exit(3);
    }
  }

  /* get the list sequence names */
  GetSeqNames(seq_list);

  /* spread the MSPs */
  while (fgets(line, 150, msp_file)) {
    if (scan_msp(&cur, line)) {
      fprintf(err, "Wrong format:\n%s", line);
      continue;
    }
    img_ct ++;
    i = GetOutputInterval(cur.query.frag.seq_name, seq_no, noof);
    if (i < 0) {
      fprintf(err, "Sequence name %s not found in the given list of sequences.  Exit.\n", cur.query.frag.seq_name);
      exit(4);
    }
    fprintf(*(output+i), "%d %d %s %d %d \n", img_ct, cur.score, cur.query.frag.seq_name, cur.query.frag.lb, cur.query.frag.rb);
    img_ct ++;
    i = GetOutputInterval(cur.sbjct.frag.seq_name, seq_no, noof);
    if (i < 0) {
      fprintf(err, "Sequence name %s not found in the given list of sequences.  Exit.\n", cur.query.frag.seq_name);
      exit(4);
    }
    fprintf(*(output+i), "%d %d %s %d %d \n", img_ct, cur.score, cur.sbjct.frag.seq_name, cur.sbjct.frag.lb, cur.sbjct.frag.rb);
  }

  fprintf(msp_no, "%d\n", img_ct/2+1);

  exit(0);
}


int GetOutputInterval(char *seq_name, int seq_no, int noof) {
  int interval;
  interval = GetSeqIndex(0, seq_no-1, seq_name)/(seq_no/noof);
  if (interval < 0) {
    printf("Sequence name %s not found in the given list of sequences.  Exit.\n", seq_name);
  }
  if (interval > noof-1) interval=noof-1; /* the last few ones */

  return interval;
}
