#ifndef __BUILD_REPEAT_FAMILIES_H__
#define __BUILD_REPEAT_FAMILIES_H__

struct posllist
{
  int this; struct posllist *next;
};
struct llist 
{ 
  int freq; int lastocc; struct llist *next; struct posllist *pos;
};

/* Note: 20 is a potential bug here.  If you choose -l 21, this code will segfault. */
struct repeatllist
{
  char value[20]; int repeatpos; struct repeatllist *next;
};

/* These are all just function forwards. */
void allocate_space();
void print_parameters();
int build_sequence(char *sequence, char *filename);
void add_rc(char *sequence);
void build_headptr(struct llist **headptr);
void trim_headptr(struct llist **headptr);
void build_all_pos(struct llist **headptr);
struct llist *find_besttmp(struct llist **headptr);
void build_pos(struct llist *besttmp);  
void extend_right();
void extend_left();
void mask_headptr(struct llist **headptr); /* removed[x]=1 means l-mer x was removed */
int maskextend_right( char mbase, int seqpos, int modelpos );
int maskextend_left( char mbase, int seqpos, int modelpos );

int hash_function(char *lmer);
int smallhash_function(char *lmer);
int lmermatch(char *lmer1, char *lmer2);
int lmermatcheither(char *lmer1, char *lmer2); /* forward or rc match */
int lmermatchrc(char *lmer1, char *lmer2); /* rc match */
int mismatches(char *lmer1, char *lmer2);
int is_degenerate(char *this);
void compute_totalbestscore_right(int y);
void print_totalbestscore_right(int y);
void compute_totalbestscore_left(int w);
void print_totalbestscore_left(int w);
void fprint_totalbestscore_left(int w);
int compute_score_right(int y, int n, int offset, char a);
int compute_score_left(int w, int n, int offset, char a);
int compute_maskscore_right(char mbase, int repeaty, int seqpos, int offset);
int compute_maskscore_left(char mbase, int repeatw, int seqpos, int offset);
int minimum_edit_distance(char *lmer, int r);
int edit_distance(char *lmer1, char *lmer2);
double compute_entropy(char *lmer);
char num_to_char(char z);
char char_to_num(char c);
char compl( char c );

#endif
