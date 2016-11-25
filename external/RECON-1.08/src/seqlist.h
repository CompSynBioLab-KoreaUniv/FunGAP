
#include "bolts.h"
#include "string.h"

#ifndef _seqlist_h
#define _seqlist_h

void GetSeqNames(FILE *);
int GetSeqIndex(int, int, char *);



void GetSeqNames(FILE *seq_list) {
  char line[256];
  int seq_ct = -1;
  char *name_start;

  /* get the number of sequences in the list, and generates an array of proper size to hold the list */
  fgets(line, 255, seq_list);
  seq_no = atoi(line);
  if (!seq_no) {
    printf("First line in the list of sequences must be an integer,\nwhich is the number of sequences in the list.  Exit\n");
    exit(2);
  }
  seq_names = (char **) malloc(seq_no*sizeof(char *));

  /* get the names of the sequences */
  while (fgets(line, 255, seq_list)) {
      seq_ct ++;
      if (seq_ct == seq_no) {
	printf("More sequence names than indicated at the begining of the file.  Exit.\n");
	exit(2);
      }
      *(seq_names+seq_ct) = (char *) malloc(NAME_LEN*sizeof(char));
      /* skip the '>' and possilbe spaces following it at the beginning of the line, in case the name is grepped directly from a fasta file */
      name_start = line;
      if (line[0] == '>') {
	name_start++;
	while (isspace(*name_start)) name_start++;
      } 
      strncpy(*(seq_names+seq_ct), name_start, NAME_LEN-1);
      /* mark the end of the string */
	 name_start = *(seq_names+seq_ct);
	 while (/**name_start != EOF &&*/ !isspace(*name_start)) name_start++;
	 *name_start = '\0';
  }

  /* we assume that the list is already sorted lexically, which is done in get-names.pl */
}




/* binary search of a sorted list of sequence names */
/* returns position in the array if found */
/* returns -1 if not found */
int GetSeqIndex(int left, int right, char *seq_name) {
  int pos, dir;

  if (left <= right) {
    pos = (left+right)/2;
    dir = strncmp(seq_name, *(seq_names+pos), NAME_LEN);
    if (dir < 0) return GetSeqIndex(left, pos-1, seq_name);
    else if (dir > 0) return GetSeqIndex(pos+1, right, seq_name);
    else return pos;
  } else {
    return -1;
  }
}


#endif
