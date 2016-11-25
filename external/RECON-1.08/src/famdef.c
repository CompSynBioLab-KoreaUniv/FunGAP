#include "ele.h"


typedef struct Dirs {
  int dir;
  struct Dirs *next;
} Dirs_t;


FAMILY_t *new_family();
void build_family(FAMILY_t *, ELE_INFO_t *, short, int *);
void recruit_partner(ELE_INFO_t *, EDGE_TREE_t *, ELE_DATA_t **, Dirs_t *, Dirs_t **);






int main (int argc, char *argv[]) {
  int i, rounds=0, ei;
  char line[35], stat;
  FAM_DATA_t *fam;
  FAMILY_t *family;
  FILE *seq_list, *ele_no, *redef_stat;
  FILE  *fam_no, *final_ele_no;
  int tot_ele;

  /* processing command line */
  if (argc == 1) {
    printf("usage: famdef seq_list\n where seq_list is the list of sequence names.\n");
    exit(1);
  }

  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  GetSeqNames(seq_list);

  ele_no = fopen("summary/redef_ele_no", "r");
  if(!ele_no) {
    printf("Can not open redef_ele_no.  Exiting\n");
    exit(2);
  }

  redef_stat = fopen("tmp/redef_stat", "r");
  if (!redef_stat) {
    printf("Can not open redef_stat.  Exiting\n");
    exit(2);
  }

  eles = fopen("summary/eles", "w");
  fams = fopen("summary/families", "w");
  fam_no = fopen("summary/fam_no", "w");
  final_ele_no = fopen("summary/ele_no", "w");

  log_file = fopen("tmp/log2", "w");

  while (fgets(line, 15, ele_no)) {
    ele_ct = atoi(line);
  }
  fclose(ele_no);

  ele_array_size = ele_ct;
  all_ele = (ELE_INFO_t **) malloc(ele_array_size*sizeof(ELE_INFO_t *));

  for (i=0; i<ele_ct; i++) {
    *(all_ele+i) = ele_info_init(i+1);
  }

  ele_ct = 0;
  while (fgets(line, 35, redef_stat)) {
    ele_ct ++;
    /*sscanf(line, "%d %c %d\n", &ei, &stat, &fu);*/
    /* for some bizzare reason, scanf doesn't work properly here. :( */
    for (i=0; i<35; i++) {
      if (line[i] == ' ') break;
    }
    ei = atoi(line);
    stat = line[i+1];
    (*(all_ele+ei-1))->stat = stat;
  }


  msp_in_mem = 0;
  msp_left = 0;
  msp_ct = 0;
  edge_ct = 0;
  edge_in_mem = 0;
  edge_left = 0;
  files_read = 0;
  clan_ct = 0;
  err_no = 0;

  /* define families */
  fam_ct = 0;
  tot_ele = 0;

  fprintf(eles, "#  fam    ele   dir  sequence    start     end\n");
  fprintf(fams, "# fam  cp-no name\n");  

  for (i=0; i<ele_ct; i++) {
    if ((*(all_ele+i))->stat != 'O' && (*(all_ele+i))->stat != 'X') {
      fprintf(log_file, "evaluating family membership of element %d\n", (*(all_ele+i))->index);
      fflush(log_file);
      if ((*(all_ele+i))->stat == 'y')	{ /* starting a new family! */
	family = new_family();
	build_family(family, *(all_ele+i), 1, &tot_ele);
      } else if ((*(all_ele+i))->stat != 'x') {
	err_no ++;
	fprintf(log_file, "ele %d %c not properly filtered\n",  (*(all_ele+i))->index, (*(all_ele+i))->stat);
	fflush(log_file);
	exit(5);
      }
    }
  }

  fprintf(fam_no, "%d\n", fam_ct);
  fclose(fam_no);

  fprintf(final_ele_no, "%d\n", tot_ele);
  fclose(final_ele_no);

  fprintf(log_file, "total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  fprintf(log_file, "%d files read, %d msps seen, %d edges seen\n", files_read, msp_ct, edge_ct);
  fprintf(log_file, "%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);


  exit(0);
}







/**************************************
 **************************************
 **                                  **
 **    CLUSTERING, YEAH BABY YEAH!   **
 **                                  **
 **************************************
 **************************************/






FAMILY_t *new_family() {
  FAMILY_t *fam;
  FAM_DATA_t *fam_data;

  fam_ct ++;
  fam = (FAMILY_t *) malloc(sizeof(FAMILY_t));
  fam->index = fam_ct;
  fam->members = NULL;
  fam->relatives = NULL;

  return fam;
}





void build_family(FAMILY_t *fam, ELE_INFO_t *ele_info, short dir, int *te) {
  /* standard breath first traverse for connected components */
  ELE_DATA_t *mem_tail, *seed, *cur_mem;
  ELE_INFO_t *cei;
  int mem_ct = 0;
  Dirs_t *dirs, *dirs_t, *cur_dir;


  ele_info->stat = 'x';
  seed = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  seed->ele_info = ele_info;
  seed->next = NULL;
  fam->members = seed;
  mem_tail = seed;
  dirs = (Dirs_t *) malloc(sizeof(Dirs_t));
  dirs->dir = 1;
  dirs->next = NULL;
  dirs_t = dirs;

  /* recruiting new members */
  cur_mem = fam->members;
  cur_dir = dirs;
  while (cur_mem) {
    mem_ct ++;
    cei = cur_mem->ele_info;
    if (cei->ele) {
      fprintf(log_file, "ele %d is in memory unexpectedly\n", cei->index);
      fflush(log_file);
      exit(3);
    } else ele_read_in(cei, 4);

    cei->to_family = fam;

    /* output the element */
    /* to the list of all elements */
    /*fprintf(eles, "%6d %6d %2d %10s %8d %8d \n", fam->index, cei->index, cei->ele->direction, cei->ele->frag.seq_name, cei->ele->frag.lb, cei->ele->frag.rb);*/
    fprintf(eles, "%6d %6d %2d %10s %8d %8d \n", fam->index, cei->index, cur_dir->dir, cei->ele->frag.seq_name, cei->ele->frag.lb, cei->ele->frag.rb);

    /* recruit its partners (that have not been recruited) */
    if (cei->ele->edges) recruit_partner(cei, cei->ele->edges, &mem_tail, cur_dir, &dirs_t);
    else if (cur_mem != seed) {
      fprintf(log_file, "error: ele %d edge tree missing\n", cei->index);
      fflush(log_file);
      exit(4);
    }
    ele_cleanup(&cei->ele);
    cur_mem = cur_mem->next;
    cur_dir = cur_dir->next;
  }

  fprintf(fams, "%d %d unknown \n", fam->index, mem_ct);
  fflush(fams);

  *te += mem_ct;

  fprintf(log_file, "new family with %d members\n", mem_ct);
  fflush(log_file);

  fam_cleanup(&fam);
}




void recruit_partner(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, ELE_DATA_t **tail, Dirs_t *cur_dir, Dirs_t **dtail){
  ELE_INFO_t *epi;
  ELE_DATA_t *new_mem;
  Dirs_t *new_dir;

  if (rt->l) recruit_partner(ele_info, rt->l, tail, cur_dir, dtail);

  if (rt->to_edge->type == 'p' || rt->to_edge->type == 'P') {
    epi = linked_ele(ele_info, rt->to_edge);
    if (epi->stat == 'y') {
      epi->stat = 'x';
      new_mem = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
      new_mem->ele_info = epi;
      new_mem->next = NULL;
      (*tail)->next = new_mem;
      *tail = new_mem;
      new_dir = (Dirs_t *) malloc(sizeof(Dirs_t));
      new_dir->dir = cur_dir->dir * rt->to_edge->direction;
      new_dir->next = NULL;
      (*dtail)->next = new_dir;
      *dtail = new_dir;
    }
  }

  if (rt->r) recruit_partner(ele_info, rt->r, tail, cur_dir, dtail);
}
