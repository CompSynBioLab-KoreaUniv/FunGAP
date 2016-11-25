#ifndef __BUILD_LMER_TABLE_H__
#define __BUILD_LMER_TABLE_H__

struct llist 
{ 
  int freq; 
  int lastplusocc; 
  int lastminusocc; 
  struct llist *next;
};

int build_sequence(char *sequence, char *filename);
void build_headptr(struct llist **headptr);
void trim_headptr(struct llist **headptr);
int count_headptr(struct llist **headptr);
void sort_headptr(struct llist **headptr, int *sortedfreq, int *sortedocc, long *sortedindex, int ngoodlmers);
void print_lmers(struct llist **headptr, int *sortedfreq, int *sortedocc, long *sortedindex, int ngoodlmers);
void q_sort(int *sortedfreq, int *sortedocc, long *sortedindex, int *order, int left, int right);

int hash_function(char *lmer);
int lmermatch(char *lmer1, char *lmer2); /* forward match */
int lmermatchrc(char *lmer1, char *lmer2); /* rc match */
int is_degenerate(char *this);
long compute_index (char *this);
void reverse_if_necessary(char *this);
char num_to_char(char z);
char char_to_num(char c);

#endif
