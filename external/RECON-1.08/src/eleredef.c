#include <fcntl.h>
#include <sys/types.h>
#include <dirent.h>

#include "ele.h"

#define CUTOFF1 0.5
#define CUTOFF2 0.9
#define MAX_IMG 1200
#define TOO_SHORT 30
#define FUDGE 2
#define MARGIN 10000
#define FLURRY 10



typedef struct img_node {
  short recorded;
  IMAGE_t *to_image;
  struct img_node *sib;
  struct img_node *children;
} IMG_NODE_t;




void report_cts();
void report_redef_stat();

ELE_DATA_t *ele_def(IMG_DATA_t **, float);
void generate_img_tree(ELEMENT_t *);
ELE_INFO_t *new_element();
void add_ele_info(ELE_INFO_t *);

void general_ele_redef(ELE_INFO_t *, IMAGE_t **);
/*void build_local_network(ELE_INFO_t *, int, ELE_DATA_t **, IMAGE_t **);
  void local_net_walk(ELE_INFO_t *, int, EDGE_TREE_t *, ELE_DATA_t **, IMAGE_t **);*/
void build_local_network(ELE_INFO_t *, ELE_DATA_t **, ELE_DATA_t **, IMAGE_t **);
void recruit(ELE_INFO_t *, EDGE_TREE_t *, ELE_DATA_t **, IMAGE_t **);
void cruise_local_net(ELE_DATA_t *, IMAGE_t **);
void local_ele_redef(ELE_INFO_t *, IMAGE_t **, int*);
void dissolve_local_network(ELE_DATA_t **);
int emptyDir(char *);
void stage_exit(int);
void dismiss_element(ELE_INFO_t *);
void remove_ele(ELE_INFO_t *);

void ele_redef(ELE_INFO_t *, IMAGE_t **);
IMG_DATA_t *img_data_sort(IMG_DATA_t *, int);
void PCP_to_TBDs(ELEMENT_t *);
BD_t *CP_cluster(CP_t *);
void CP_sort(CP_t **);
void BD_sort(BD_t **);
int span(ELEMENT_t *, int32_t);
void TBD_merge(ELEMENT_t *);

void dissect(ELE_INFO_t *);
int too_short(FRAG_t *);
MSP_t *add_msp(MSP_t *);
void register_image(IMAGE_t *, ELEMENT_t *);
void put_image(IMG_DATA_t **, IMAGE_t *);
void dump_image(IMAGE_t *);
void remove_image(IMAGE_t *);
void combo_update(ELE_INFO_t *);
void combo_edge_update(ELE_INFO_t *, EDGE_TREE_t **);
void CP_clean(CP_t **, ELE_INFO_t *);

void edges_and_cps(ELE_INFO_t *, IMAGE_t **);
int full_length(IMAGE_t *, float);
void add_CP(CP_t **, int32_t, ELE_INFO_t *);
void add_edge(ELE_INFO_t *, ELE_INFO_t *, char, int32_t, short);
void adjust_edge_tree(ELE_INFO_t *);
int charge_edge_array(EDGE_t **, EDGE_TREE_t *, int);
int consis_tree_build(IMG_NODE_t *, IMAGE_t *, int);
int consis(IMAGE_t *, IMAGE_t *, float);
IMG_NODE_t **node_entry(IMG_NODE_t **);
void consis_tree_free(IMG_NODE_t *);
/*int find_prim(struct img_node *, float, int32_t, int32_t, short, int32_t *, short *);*/
int find_prim(IMG_NODE_t *, float, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int *, int32_t *, short *);

void combo_output(ELE_INFO_t *);
void obs_output(ELE_INFO_t *);
void fprint_ele_obs(FILE *, ELE_INFO_t *);

int frag_cmp(const void *, const void *);
int partner_cmp(const void *, const void *);
int CP_cmp(const void *, const void *);
int BD_cmp(const void *, const void *);
int fam_cmp(const void *, const void *);



int main (int argc, char *argv[]) {
  ELE_INFO_t *cur_ele_info;
  int i, ele_march, ei, rounds=0, start;
  char line[35], stat;
  short fu, to_march;
  IMAGE_t **img_ptr;
  FILE *ele_no, *msp_no, *edge_no, *size_list, *new_stat;
  FILE *seq_list;

  /* processing command line */
  if (argc == 1) {
    printf("usage: redef seq_list start clan_ct\n where seq_list is the list of sequence names, start is the index of the element to start redefining, and clan_ct is the number to start counting the number of clans.  The latter two are optional.\n");
    exit(1);
  }

  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  GetSeqNames(seq_list);

  if (argc > 2) start = atoi(argv[2]) - 1;
  else start = 0;
  if (argc > 3) clan_ct = atoi(argv[3]);
  else clan_ct = 0;

  ele_no = fopen("summary/naive_ele_no", "r");
  if(!ele_no) {
    printf("Can not open naive_ele_no.  Exiting\n");
    exit(1);
  }
  msp_no = fopen("summary/redef_msp_no", "r");
  if (!msp_no) msp_no = fopen("summary/ori_msp_no", "r");
  if (!msp_no) {
    printf("Can not open msp_no.  Exiting\n");
    exit(1);
  }

  // May or may not exist.  All uses of edge_no first check for
  // its existence.
  edge_no = fopen("summary/naive_edge_no", "r");

  // May or may not exist.  All uses are checked.
  size_list = fopen("tmp/size_list", "r");

  // May or may not exist.  All uses are checked.
  new_stat = fopen("tmp2/redef_stat", "r");

  new_msps = fopen("summary/new_msps", "a");
  if ( !new_msps )
  {
    printf("Can not open summary/new_msps for writing! Exiting\n");
    exit(1);
  }
  unproc = fopen("summary/unproc", "a");
  if ( !unproc )
  {
    printf("Can not open summary/unproc for writing! Exiting\n");
    exit(1);
  }
  combo = fopen("summary/combo", "a");
  if ( !combo )
  {
    printf("Can not open summary/combo for writing! Exiting\n");
    exit(1);
  }
  obs = fopen("summary/obsolete", "a");
  if ( !obs )
  {
    printf("Can not open summary/obsolete for writing! Exiting\n");
    exit(1);
  }
  log_file = fopen("tmp2/log", "a");
  if ( !log_file )
  {
    printf("Can not open tmp2/log for writing! Exiting\n");
    exit(1);
  }

  while (fgets(line, 15, ele_no)) {
    ele_ct = atoi(line);
  }
  fclose(ele_no);

  while (fgets(line, 15, msp_no)) {
    msp_index = atoi(line) - 1;
  }
  fclose(msp_no);

  if (edge_no) {
    while (fgets(line, 15, edge_no)) {
      edge_index = atoi(line) - 1;
    }
  } else edge_index = -1;

  ele_array_size = 2*ele_ct;
  all_ele = (ELE_INFO_t **) malloc(ele_array_size*sizeof(ELE_INFO_t *));
  for (i=0; i<ele_array_size; i++) {
    *(all_ele+i) = ele_info_init(i+1);
  }
  ele_info_data = NULL;

  /* get rid of large tandem reps */
  if (size_list && !new_stat) outthrow_big_tandems(size_list);
  fclose(unproc);

  /* set file_updated / stat for ele_info's */
  if (new_stat) {
    ele_ct = 0;
    while (fgets(line, 35, new_stat)) {
      ele_ct ++;
      /*sscanf(line, "%d %c %d\n", &ei, &stat, &fu);*/
      /* for some bizzare reason, scanf doesn't work properly here. :( */
      for (i=0; i<35; i++) {
	if (line[i] == ' ') break;
      }
      ei = atoi(line);
      stat = line[i+1];
      fu = atoi(&line[i+3]);
      if (ei<=ele_array_size) {
	(*(all_ele+ei-1))->stat = stat;
	(*(all_ele+ei-1))->file_updated = fu;
      } else {
	cur_ele_info = ele_info_init(ei);
	cur_ele_info->stat = stat;
	cur_ele_info->file_updated = fu;
	add_ele_info(cur_ele_info);
      }
    }
  }

  msp_in_mem = 0;
  msp_left = 0;
  msp_ct = 0;
  edge_ct = 0;
  edge_in_mem = 0;
  edge_left = 0;
  files_read = 0;
  err_no = 0;

  /* re-define elements using the syntopy algorithm, and build edges */
  img_ptr = (IMAGE_t **) malloc(MAX_IMG*sizeof(IMAGE_t *));

  if ( ! img_ptr )
  {
    printf("eleredef: Error! Could not allocate memory for img_ptr: %d bytes requested\n", ( MAX_IMG*sizeof(IMAGE_t *) ) );
    exit(-1);
  }

  to_march = 1;
  while (to_march) {
    to_march = 0;
    rounds ++;
    for (i=start; i<ele_ct && i<ele_array_size; i++) {
      fprintf(log_file, "evaluating definition of element %d\n", (*(all_ele+i))->index);
      fflush(log_file);
      if ((*(all_ele+i))->stat == 'z' || (*(all_ele+i))->stat == 'w' || (*(all_ele+i))->stat == 't') {
	to_march = 1;
	general_ele_redef(*(all_ele+i), img_ptr);
#if 0
	report_redef_stat();
#endif
      } /*else if ((*(all_ele+i))->stat == 'O' && !(*(all_ele+i))->file_updated) spit_out_ele(*(all_ele+i));*/
    }

    start = 0;
    
    cur_ele_info = ele_info_data;
    while(cur_ele_info) {
      fprintf(log_file, "evaluating definition of element %d\n", cur_ele_info->index);
      fflush(log_file);
      if (cur_ele_info->stat == 'z' || cur_ele_info->stat == 'w' || cur_ele_info->stat == 't') {
	to_march = 1;
	general_ele_redef(cur_ele_info, img_ptr);
#if 0
	report_redef_stat();
#endif
      } /*else if (cur_ele_info->stat == 'O' && !cur_ele_info->file_updated) spit_out_ele(cur_ele_info);*/
      cur_ele_info = cur_ele_info->next;
    }
  }

  report_cts();
  report_redef_stat();
  free(img_ptr);

#if 0
  for (i=start; i<ele_ct && i<ele_array_size; i++) {
    if ((*(all_ele+i))->stat == 'X') remove_ele(*(all_ele+i));
  }
    
  cur_ele_info = ele_info_data;
  while(cur_ele_info) {
    if (cur_ele_info->stat == 'X') remove_ele(cur_ele_info);
    cur_ele_info = cur_ele_info->next;
  }
#endif

  fprintf(log_file, "total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  fprintf(log_file, "%d rounds, %d files read, %d msps seen, %d edges seen\n", rounds, files_read, msp_ct, edge_ct);
  fprintf(log_file, "%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);

  exit(0);
}



void report_cts() {
  FILE *fp;

  fp = fopen("summary/redef_ele_no", "w");
  fprintf(fp, "%d\n", ele_ct);
  fclose(fp);
}




void report_redef_stat() {
  int i, ele_march = ele_ct < ele_array_size ? ele_ct : ele_array_size;
  FILE *redef_stat, *fp;
  ELE_INFO_t *cur_ele_info = ele_info_data;

  redef_stat = fopen("tmp2/redef_stat", "w");
  for (i=0; i<ele_march; i++) {
    fprintf(redef_stat, "%d %c %d\n", (*(all_ele+i))->index, (*(all_ele+i))->stat, (*(all_ele+i))->file_updated);
  }
  while(cur_ele_info) {
    fprintf(redef_stat, "%d %c %d\n", cur_ele_info->index, cur_ele_info->stat, cur_ele_info->file_updated);
    cur_ele_info = cur_ele_info->next;
  }
  fclose(redef_stat);

  fp = fopen("summary/redef_msp_no", "w");
  fprintf(fp, "%d\n", msp_index+1);
  fclose(fp);

  fp = fopen("summary/naive_edge_no", "w");
  fprintf(fp, "%d\n", edge_index+1);
  fclose(fp);
}


/***********************************
 ***********************************
 **                               **
 **       DEFINING  ELEMENTS      **
 **                               **
 ***********************************
 ***********************************/


ELE_DATA_t *ele_def(IMG_DATA_t **img_data_p, float cutoff) {
  /* this function assumes that the images are SORTED according to frag_cmp.
     defining elements is in essence partitioning images in the original list.
     cur_img_data marks the current image under inspection;
     *img_data_p and prev_img_data are used to maintain the un-participated list;
     *img_data_p marks the first (un-participated) image in the remaining list;
     prev_img_data marks the last un-participated image proceding cur_img_data, so that when cur_img_data is removed from the origial list and participated into an element, we know where the tail of the unparticipated list, used to be pointed by cur_img_data->next, should be appended to */


  IMG_DATA_t *cur_img_data, *prev_img_data;
  ELEMENT_t *ele_tmp;
  ELE_INFO_t *ele_info_tmp;
  ELE_DATA_t *ele_data=NULL, *ele_data_tmp;
  short ritetime;
  int ct;

 while (*img_data_p) {

  cur_img_data = *img_data_p;
  prev_img_data = NULL;
  ritetime = 1; /* ritetime marks the right time to start a new element */

  while(cur_img_data != NULL) {
    if (ritetime) { /* starting a new element */
      ele_info_tmp = new_element();
      ele_tmp = ele_info_tmp->ele;
      /* start the definition of the element */
      ele_tmp->img_no = 1;
      ele_tmp->frag = cur_img_data->to_image->frag; /* lst and rst included */
      /* kick cur_img_data out of the original list */
      /* notice that whenever we're here, cur_img_data points to the first image in the remaining list, which is also pointed to by *img_data_p */
      *img_data_p = cur_img_data->next;
      /* updat the current image and move it to the image list of ele_tmp */
      /* notice that here it's always the first one for ele_tmp */
      cur_img_data->to_image->ele_info = ele_info_tmp;
      cur_img_data->next = NULL;
      ele_tmp->to_img_data = cur_img_data;
      /* move pointer cur_img_data to the next one */
      cur_img_data = *img_data_p;
      ritetime = 0;
      continue;
    }
    if (/*!strncmp(ele_tmp->frag.seq_name, cur_img_data->to_image->frag.seq_name, NAME_LEN)*/ ele_tmp->frag.seq_name == cur_img_data->to_image->frag.seq_name && ele_tmp->frag.rb - cur_img_data->to_image->frag.lb > 10) { /*checking possible images */
      if (sing_cov(&ele_tmp->frag, &cur_img_data->to_image->frag, cutoff)) { /* a good image */
	ele_tmp->img_no ++;
	/* move cur_img_data from the original list */
	/* notice that it's a bit more complicated than above */
	if (prev_img_data != NULL) {
	  prev_img_data->next = cur_img_data->next;
	} else {
	  *img_data_p = cur_img_data->next;
	}
	/* update the current image and move it to the image_list of ele_tmp */
	cur_img_data->to_image->ele_info = ele_info_tmp;
	cur_img_data->next = ele_tmp->to_img_data;
	ele_tmp->to_img_data = cur_img_data;

	/* updating definition of the element if necessary, and move cur_img_data to the right place */
	/* notice that it is a bit complicated than above */
	if (ele_tmp->frag.rb < cur_img_data->to_image->frag.rb) { 
	  ele_tmp->frag.rb = cur_img_data->to_image->frag.rb;
	  /* ele_tmp->frag.rst = cur_img_data->to_image->frag.rst; */
	  /* time to go back and check if previously unqualified images are now good to be participated to the updated element */
	  cur_img_data = *img_data_p;
	  prev_img_data = NULL;
	} else {
	  if (prev_img_data != NULL) cur_img_data = prev_img_data->next;
	  else cur_img_data = *img_data_p;
	}
      } else { /* current image not good for the ele_tmp, keep it in the remaining list and move on */
	prev_img_data = cur_img_data;
	cur_img_data = cur_img_data->next;
      }
    } else { /* time to finish defining of ele_tmp */
      /* add the element to the list of elements */
      ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
      ele_data_tmp->ele_info = ele_info_tmp;
      ele_data_tmp->next = ele_data;
      ele_data = ele_data_tmp;
      /*      if (!strncmp(ele_tmp->frag.seq_name, target, NAME_LEN)) printf("%d %s %d %d\n", ele_tmp->index, ele_tmp->frag.seq_name, ele_tmp->frag.lb, ele_tmp->frag.rb); */
      /* time to go back to the head of the remaining list and start a new element */
      ritetime = 1;
      cur_img_data = *img_data_p;
      prev_img_data = NULL;
    }
  }
  /* add the last element to the list of elements */
  ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  ele_data_tmp->ele_info = ele_info_tmp;
  ele_data_tmp->next = ele_data;
  ele_data = ele_data_tmp;

 }

  ele_data_tmp = ele_data;
  while (ele_data_tmp) {
    if (ele_data_tmp->ele_info->ele->frag.lb > ele_data_tmp->ele_info->ele->frag.rb) {
      fprintf(log_file, "error:  ele %d reversed after ele_def\n", ele_data_tmp->ele_info->index);
      fflush(log_file);
      exit(3);
    }
    generate_img_tree(ele_data_tmp->ele_info->ele);
    /*ct = count_img_nodes(ele_data_tmp->ele_info->ele->to_img_tree);
    if (ct != ele_data_tmp->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  trouble generating image tree for redefed offsprings: ele %d %d %d\n", ele_data_tmp->ele_info->index, ele_data_tmp->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
    }*/
    ele_data_tmp = ele_data_tmp->next;
  }

  return ele_data;
}




void generate_img_tree(ELEMENT_t *ele) {
  IMAGE_t **img_array = (IMAGE_t **) malloc(ele->img_no*sizeof(IMAGE_t *));
  IMG_DATA_t *cur;
  int ct=0;

  cur = ele->to_img_data;
  while (cur) {
    *(img_array+ct) = cur->to_image;
    ct ++;
    cur = cur->next;
  }

  qsort(img_array, ct, sizeof(IMAGE_t *), img_index_cmp);

  build_img_tree(&ele->to_img_tree, img_array, 0, ct-1);

  free(img_array);
}





ELE_INFO_t *new_element() {
  ELE_INFO_t *ele_info_tmp;

  ele_ct ++;

  if (ele_ct <= ele_array_size) ele_info_tmp = *(all_ele+ele_ct-1);
  else {
    ele_info_tmp = ele_info_init(ele_ct);
    add_ele_info(ele_info_tmp);
  }

  ele_info_tmp->ele = ele_init(ele_ct);
  ele_info_tmp->file_updated = 1;

  return ele_info_tmp;
}




void add_ele_info(ELE_INFO_t *ele_info_tmp) {

  if (!ele_info_data) ele_info_data = ele_info_tmp;
  else ele_info_tail->next = ele_info_tmp;

  ele_info_tail = ele_info_tmp;
}





/***********************************
 ***********************************
 **                               **
 **    REDEFINING THE ELEMENTS    **
 **                               **
 ***********************************
 ***********************************/




/**********************************************
 * Organizing the traverse of the local graph *
 **********************************************/





void general_ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
  ELE_DATA_t *local_net=NULL, *local_net_tail=NULL;

  if (!ele_info->ele) ele_read_in(ele_info, 1);

  /* an element could have lost all its images during the re-evaluation of its partners */
  if (!ele_info->ele->img_no) {
    combo_update(ele_info);
    remove_ele(ele_info);
    ele_cleanup(&ele_info->ele);
  } else {
    /* redefining cur_ele and its adjacent neighbors */
    fprintf(new_msps, "new clan for ele %d\n", ele_info->index);
    fflush(new_msps);
    fprintf(combo, "new clan for ele %d\n", ele_info->index);
    fflush(combo);
    fprintf(obs, "new clan for ele %d\n", ele_info->index);
    fflush(obs);
    /* set up the local network */
    clan_ct ++;
    clan_size = 0;
    clan_core_size = 0;
    fprintf(log_file, "new clan: %d for ele %d\n", clan_ct, ele_info->index);
    fflush(log_file);
    build_local_network(ele_info, &local_net, &local_net_tail, img_ptr);
    if ( ele_info->index == 362 || ele_info->index == 10 )
    {
      printf("Printing local net for 362\n");
      print_ele_data( local_net );
    }
    fprintf(log_file, "clan size: %d, clan core size: %d\n", clan_size, clan_core_size);
    fflush(log_file);
    /* redefining elements in the local network */
    /* local_ele_redef(ele_info, img_ptr); */
    cruise_local_net(local_net, img_ptr);
    /* clearing up the local network */
    fflush(new_msps);
    fflush(combo);
    fflush(obs);
    dissolve_local_network(&local_net);
  }
}




#if 0
void build_local_network(ELE_INFO_t *ele_info, int level, ELE_DATA_t **net_p, IMAGE_t **img_ptr) {
  ELEMENT_t *ele;
  ELE_DATA_t *que;

  clan_size ++;
  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = *net_p;
  *net_p = que;

  if (!ele_info->ele) ele_read_in(ele_info, 1);
  ele = ele_info->ele;

  if (level < DEPTH) {
    ele->l_hold = 1;
    clan_core_size ++;
    /*fprintf(log_file, "    %d\n", ele_info->index);
      fflush(log_file);*/
    if (ele_info->stat != 'v') {
      if (ele_info->stat == 'z') edges_and_cps(ele_info, img_ptr);
      level ++;
      if (ele->edges) local_net_walk(ele_info, level, ele->edges, net_p, img_ptr);
    } else {
      if (ele->edges) local_net_walk(ele_info, DEPTH, ele->edges, net_p, img_ptr);
    }
  } else ele->l_hold = 2;
}




void local_net_walk(ELE_INFO_t *ele_info, int level, EDGE_TREE_t *edge_node, ELE_DATA_t **net_p, IMAGE_t **img_ptr) {
  ELE_INFO_t *epi;

  if (edge_node->l) local_net_walk(ele_info, level, edge_node->l, net_p, img_ptr);

  /*if (edge_node->to_edge->type == 'p' || edge_node->to_edge->type == 's') {*/
    epi = linked_ele(ele_info, edge_node->to_edge);
    if (!epi->ele) ele_read_in(epi, 1);
    if (!epi->ele->l_hold) build_local_network(epi, level, net_p, img_ptr);
    else if (epi->ele->l_hold == 2 && level < DEPTH) {
      epi->ele->l_hold = 1;
      clan_core_size ++;
      /*fprintf(log_file, "    %d\n", epi->index);
	fflush(log_file);*/
      if (epi->stat == 'z') edges_and_cps(epi, img_ptr);
      if (epi->ele->edges) local_net_walk(epi, level+1, epi->ele->edges, net_p, img_ptr);
    }
  /*}*/
  
  if (edge_node->r) local_net_walk(ele_info, level, edge_node->r, net_p, img_ptr);
}
#endif





void build_local_network(ELE_INFO_t *ele_info, ELE_DATA_t **net_p, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELEMENT_t *ele;
  ELE_DATA_t *que;

  /* seed the network with the first element */
  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = NULL;
  *net_p = que;
  *net_tail_p = que;

  ele_info->ele->l_hold = 1;

  /* breadth first search */
  while (que) {
    clan_size ++;
    if (que->ele_info->ele->l_hold <= DEPTH) {
      clan_core_size ++;
      if (que->ele_info->stat == 'z') edges_and_cps(que->ele_info, img_ptr);
      if (que->ele_info->ele->edges) recruit(que->ele_info, que->ele_info->ele->edges, net_tail_p, img_ptr);
      /* else {
	fprintf(log_file, "error: ele %d edge tree missing.  Exit.\n", que->ele_info->index);
	fflush(log_file);
	exit(4);
      }*/
    }
    que = que->next;
  }
}




void recruit(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELE_INFO_t *epi;
  ELE_DATA_t *member;

  if (rt->l) recruit(ele_info, rt->l, net_tail_p, img_ptr);

  epi = linked_ele(ele_info, rt->to_edge);
  if (!epi->ele) ele_read_in(epi, 1);
  if (!epi->ele->l_hold) {
    epi->ele->l_hold = ele_info->ele->l_hold + 1;
    if (epi->stat == 'v' && epi->ele->l_hold < DEPTH) epi->ele->l_hold = DEPTH;
    member = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
    member->ele_info = epi;
    member->next = NULL;
    (*net_tail_p)->next = member;
    *net_tail_p = member;
  }

  if (rt->r) recruit(ele_info, rt->r, net_tail_p, img_ptr);
}




void cruise_local_net(ELE_DATA_t *local_net, IMAGE_t **img_ptr) {
  ELE_DATA_t *que;
  int to_march = 1;

  while (to_march) {
    to_march = 0;
    que = local_net;
    while (que) {
	if (que->ele_info->ele->l_hold <= DEPTH) {
	  local_ele_redef(que->ele_info, img_ptr, &to_march);
	}
	que = que->next;
    }
  }
}




void local_ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr, int *march_p) {
    ELE_DATA_t *new_ele_data;
    ELEMENT_t *ele = ele_info->ele;
    
    if (ele->redef != NULL) {
      new_ele_data = ele->redef;
      while (new_ele_data != NULL) {
        local_ele_redef(new_ele_data->ele_info, img_ptr, march_p);
        new_ele_data = new_ele_data->next;
      }
    } else {
	if (ele_info->stat != 'v' && ele_info->stat != 'X') {
	  if (ele->PCP) {
  	    *march_p = 1;
	    ele_redef(ele_info, img_ptr);
	  } else ele_info->stat = 'v';
	}
    }
}   





void dissolve_local_network(ELE_DATA_t **net_p) {
  ELE_DATA_t *que, *new_ele_data;
  int i, ele_left=0;
  char *command;
  /*MSP_DATA_t *md;*/

  int j, in_clan;
  FILE *fp;
  ELE_INFO_t *cur_ele_info;

  /* output results to tmp2/clan, which holds results temporily and clean up memory */
  que = *net_p;
  while (que) {
    dismiss_element(que->ele_info);
    que = que->next;
  }

  if (msp_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d msps total, %d seen, %d left in memory\n", msp_index+1, msp_ct, msp_in_mem);
    fflush(log_file);
    msp_left += msp_in_mem;
    /*md = all_msps;
    while (md) {
      if (md->to_msp) {
	fprint_msp(log_file, md->to_msp);
	fflush(log_file);
	MSP_free(md->to_msp);
      }
      md = md->next;
    }
    if (msp_in_mem) {
      fprintf(log_file, "error: further discrepency: %d msps left\n", msp_in_mem);
      msp_in_mem = 0;
    }*/
  }

  if (edge_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d edges total, %d seen, %d left in memory\n", edge_index+1, edge_ct, edge_in_mem);
    fflush(log_file);
    edge_left += edge_in_mem;
    edge_in_mem = 0;
  }

#if 0
  for (i=0; i<ele_ct; i++) {
    if ((*(all_ele+i))->ele) {
      ele_left ++;
      err_no ++;
      fprintf(log_file, "error:  element %d not cleaned from memory. stat = %c\n", i+1, (*(all_ele+i))->stat);
      fflush(log_file);
      if ((*(all_ele+i))->stat != 'O' && (*(all_ele+i))->stat != 'X') ele_write_out(*(all_ele+i), 1);
      ele_cleanup(&(*(all_ele+i))->ele);
    }
  }
  if (ele_left) {
    err_no ++;
    fprintf(log_file, "error:  %d elements left in memory\n", ele_left);
    fflush(log_file);
    /*system(command);*/
  }
#endif

  if (err_no) exit(3);

#if 0
  /* move files from tmp2/clan, i.e. finalize the result */
  command  = (char *) malloc(50*sizeof(char));

  sprintf(command, "tmp2/clan/combos");
  if (!emptyDir(command)) {
    sprintf(command, "mv -f tmp2/clan/combos combos/clan%d", clan_ct);
    if (system(command)) {
      fprintf(log_file, "error moving tmp2/clan/combos/\n");
      fflush(log_file);
      stage_exit(6);
    } else {
      sprintf(command, "mkdir tmp2/clan/combos");
      if (system(command)) {
	fprintf(log_file, "error generating dir tmp2/clan/combos/\n");
	fflush(log_file);
	stage_exit(6);
      }
    }
  }

  sprintf(command, "tmp2/clan/obsolete");
  if (!emptyDir(command)) {
    sprintf(command, "mv -f tmp2/clan/obsolete obsolete/clan%d", clan_ct);
    if (system(command)) {
      fprintf(log_file, "error moving tmp2/clan/obsolete/\n");
      fflush(log_file);
      stage_exit(6);
    } else {
      sprintf(command, "mkdir tmp2/clan/obsolete");
      if (system(command)) {
	fprintf(log_file, "error generating dir tmp2/clan/obsolete\n");
	fflush(log_file);
	stage_exit(6);
      }
    }
  }

  sprintf(command, "mv -f tmp2/clan/e* tmp2/.");
  if (system(command)) {
    /*fprintf(log_file, "error moving files from tmp2/clan/\n");
    fflush(log_file);
    stage_exit(6);*/
    j = ele_ct < ele_array_size ? ele_ct : ele_array_size;
    for (i=0; i<=j; i++) {
      sprintf(command, "tmp2/clan/e%d", i+1);
      if (fp = fopen(command, "r")) in_clan = 1;
      else in_clan = 0;
      fclose(fp);
      if (in_clan) {
	sprintf(command, "mv -f tmp2/clan/e%d tmp2/.", i+1);
	if (system(command)) {
	  fprintf(log_file, "error doing %s\n", command);
	  fflush(log_file);
	  stage_exit(6);
	}
      }
    }
    cur_ele_info = ele_info_data;
    while (cur_ele_info) {
      sprintf(command, "tmp2/clan/e%d", cur_ele_info->index);
      if (fp = fopen(command, "r")) in_clan = 1;
      else in_clan = 0;
      fclose(fp);
      if (in_clan) {
	sprintf(command, "mv -f tmp2/clan/e%d tmp2/.", cur_ele_info->index);
	if (system(command)) {
	  fprintf(log_file, "error doing %s\n", command);
	  fflush(log_file);
	  stage_exit(6);
	}
      }
      cur_ele_info = cur_ele_info->next;
    }
  }

  free(command);  
#endif

  /* remove arkaine ele files, now done in main()
  que = *net_p;
  while (que) {
    if (que->ele_info->stat == 'X') remove_ele(que->ele_info);
    que = que->next;
  } */

  ele_data_free(net_p);
}



#if 0
int emptyDir(char *dir_name) {
  int res = 1;
  DIR *cur_dir;
  struct dirent *dp;

  cur_dir = opendir(dir_name);
  while (cur_dir && res && (dp = readdir(cur_dir))) {
    if (strcmp(dp->d_name, ".") && strcmp(dp->d_name, "..")) res = 0;
  }
  closedir(cur_dir);
  return res;
}




void stage_exit(int code) {
  char *command = (char *) malloc(50*sizeof(char));

  sprintf(command, "mv -f tmp2/redef_stat tmp2/redef_stat_prev");
  if (system(command)) {
    fprintf(log_file, "error in moving tmp2/redef_stat\n");
    fflush(log_file);
    exit(code);
  }

  sprintf(command, "mv -f tmp2/msp_no tmp2/msp_no_prev");
  if (system(command)) {
    fprintf(log_file, "error in moving tmp2/msp_no\n");
    fflush(log_file);
    exit(code);
  }

  sprintf(command, "mv -f tmp2/edge_no tmp2/edge_no_prev");
  if (system(command)) {
    fprintf(log_file, "error in moving tmp2/edge_no\n");
    fflush(log_file);
    exit(code);
  }

  report_redef_stat();
  exit(code);
}
#endif



void dismiss_element(ELE_INFO_t *ele_info) {
    ELE_DATA_t *new_ele_data;

    if (ele_info->ele->redef) {
      new_ele_data = ele_info->ele->redef;
      while (new_ele_data) {
        dismiss_element(new_ele_data->ele_info);
	new_ele_data = new_ele_data->next;
      }
    } 
    /*printf("clearing ele %d\n", ele_info->index);*/
    ele_info->ele->l_hold = 0;
    if (ele_info->stat != 'X') {
      if (!ele_info->ele->img_no) {
	ele_info->stat = 'X';
	if (ele_info->ele->redef) combo_output(ele_info);
	else obs_output(ele_info);
      } else ele_write_out(ele_info, 1);
    }
    ele_cleanup(&ele_info->ele);
}




void remove_ele(ELE_INFO_t *ele_info) {
    char *command;

    command = (char *) malloc(50*sizeof(char));
    sprintf(command, "rm -f tmp2/e%d", ele_info->index);
    if (system(command)) {
      fprintf(log_file, "error removing file tmp2/e%d\n", ele_info->index);
      fflush(log_file);
      exit(6);
    }
    free(command);
}






/******************************
 * Redefining a given element *
 ******************************/



void ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
    ELE_DATA_t *new_ele_data;
    ELEMENT_t *cur_ele = ele_info->ele;
    BD_t *pbd;
    CP_t *cur;
    short to_dissect=0;
    IMG_DATA_t *cur_img_data;

    if (!cur_ele->to_img_data) listify(cur_ele->to_img_tree, &cur_ele->to_img_data);
    cur_ele->to_img_data = img_data_sort(cur_ele->to_img_data, cur_ele->img_no);

    /* find TBDs for the element */
    if (cur_ele->PCP) {
      PCP_to_TBDs(cur_ele);
    }
    if (cur_ele->TBD) {
      TBD_merge(cur_ele);
    }
    if (cur_ele->TBD) {
#if 0
      to_dissect = 1;
      if (cur_ele->TBD->bd - cur_ele->frag.lb <= FLURRY || cur_ele->TBD->bd - cur_ele->frag.rb >= -FLURRY) {
	if (!cur_ele->TBD->next) {
	  to_dissect = 0;
	  if (cur_ele->TBD->bd - cur_ele->frag.lb <= FLURRY) cur_ele->frag.lb = cur_ele->TBD->bd;
	  else cur_ele->frag.rb = cur_ele->TBD->bd;
	  BD_free(&cur_ele->TBD);
	} else {
	  if (cur_ele->TBD->next->bd - cur_ele->frag.rb >= -FLURRY) {
	    if (!cur_ele->TBD->next->next) {
	      to_dissect = 0;
	      cur_ele->frag.lb = cur_ele->TBD->bd;
	      cur_ele->frag.rb = cur_ele->TBD->next->bd;
	      BD_free(&cur_ele->TBD);
	    }
	  }
	}
      }
#endif
      pbd = cur_ele->TBD;
      while (pbd) {
	/*if (pbd->bd-cur_ele->frag.lb<-(FLURRY+5) || pbd->bd-cur_ele->frag.rb>(FLURRY+5)) {
	  fprintf(log_file, "error:  ele %d TBD outside ele bound\n", ele_info->index);
	  fflush(log_file);
	  exit(3);
	}*/
	if (pbd->bd-cur_ele->frag.lb>FLURRY && pbd->bd-cur_ele->frag.rb<-FLURRY) {
	  to_dissect ++;
	}
	pbd = pbd->next;
      }

      if (!to_dissect) {
	pbd = cur_ele->TBD;
	while (pbd) {
	  if (/*pbd->bd-cur_ele->frag.lb>0 &&*/ pbd->bd-cur_ele->frag.lb<=FLURRY) {
	    cur_ele->frag.lb = pbd->bd;
	  } else if (/*pbd->bd-cur_ele->frag.rb<0 &&*/ pbd->bd-cur_ele->frag.rb>=-FLURRY) {
	    cur_ele->frag.rb = pbd->bd;
	  }
	  pbd = pbd->next;
	}
	BD_free(&cur_ele->TBD);
      } else {
	/* dissect all images according to the TBDs */
	dissect(ele_info);
	/* redefine elements according to the dissected images, if any left */
	if (cur_ele->img_no) {
	  cur_ele->to_img_data = img_data_sort(cur_ele->to_img_data, cur_ele->img_no);
	  cur_ele->redef = ele_def(&cur_ele->to_img_data, CUTOFF1);
	}
	/* clear unnecessary memory, 'v'->'w' and update CPs for w's */
	combo_update(ele_info);
	/* ele_redef() for offspring elements.  we finish redefinition of all offsprings before pulling in partners, because the combo may have > 1 copy of the same family, which means other members' CPs won't be fully updated until all offsprings are processed */
	new_ele_data = cur_ele->redef;
	while (new_ele_data != NULL) {
	  new_ele_data->ele_info->ele->update = 1;
	  new_ele_data->ele_info->ele->l_hold = cur_ele->l_hold;
	  if (new_ele_data->ele_info->ele->to_img_tree) {
	    edges_and_cps(new_ele_data->ele_info, img_ptr);
	    if (new_ele_data->ele_info->ele->PCP) ele_redef(new_ele_data->ele_info, img_ptr);
	    else new_ele_data->ele_info->stat = 'v';
	  } else {
	   if (new_ele_data->ele_info->ele->img_no) {
	    err_no ++;
	    fprintf(log_file, "error:  image tree missing in newly redefed offspring %d\n", new_ele_data->ele_info->index);
	    fflush(log_file);
	    exit(2);
	   } else {
	    combo_output(new_ele_data->ele_info);
	   }
	  }
	  new_ele_data = new_ele_data->next;
	}
      }
    }
    if (cur_ele->img_no > 0) {
      ele_info->stat = 'v';
    } else if (!to_dissect) combo_update(ele_info);

    /*img_data_free(&cur_ele->to_img_data);*/
}





IMG_DATA_t *img_data_sort(IMG_DATA_t *img_data, int ct) {
  IMG_DATA_t **img_data_ptr;
  int i;

  if (!ct) return NULL;
  img_data_ptr = (IMG_DATA_t **) malloc(ct*sizeof(IMG_DATA_t *));
  for (i=0; i<ct; i++) {
    *(img_data_ptr+i) = img_data;
    img_data = img_data->next;
  }
  qsort(img_data_ptr, ct, sizeof(IMAGE_t *), frag_cmp);
  img_data = *img_data_ptr;
  ct --;
  for (i=0; i<ct; i++) {
    (*(img_data_ptr+i))->next = *(img_data_ptr+i+1);
  }
  (*(img_data_ptr+ct))->next = NULL;
  free(img_data_ptr);
  return img_data;
}




/***************************
 * Determine possible TBDs *
 ***************************/



void PCP_to_TBDs(ELEMENT_t *ele) {
  int s = 0, left;
  BD_t *pbd_tmp, *pbd_prev, *pbd, *pbds;
  CP_t *cp;

  /* sort the PCP list according to cp */
  CP_sort(&ele->PCP);
  /*clustering the PCPs into PBDs */
  pbds = CP_cluster(ele->PCP);
  /* identify TBP from PBDs */
  /* TBDs are removed from PBDs, what is left in PBDs are those unsuccessful ones */
  pbd_tmp = pbds;
  pbd_prev = NULL;
  while (pbd_tmp != NULL) {
    /* s is the KEY! */
    s = span(ele, pbd_tmp->bd);
    if (pbd_tmp->support >= s) {
	if (pbd_prev == NULL) {
	  pbds = pbd_tmp->next;
	  pbd_tmp->next = ele->TBD;
	  ele->TBD = pbd_tmp;
	  pbd_tmp = pbds;
	} else {
	  pbd_prev->next = pbd_tmp->next;
	  pbd_tmp->next = ele->TBD;
	  ele->TBD = pbd_tmp;
	  pbd_tmp = pbd_prev->next;
	}
    } else {
      pbd_prev = pbd_tmp;
      pbd_tmp = pbd_tmp->next;
    }
  }

  BD_free(&pbds);
}



#if 0
BD_t *CP_cluster(CP_t *cps) {
  int32_t first = cps->cp, last = cps->cp, sum = 0;
  int cpct = 0;
  BD_t *bds = NULL, *bd_tmp;

  while (cps != NULL) {
    if (cps->cp - last <= 10) {
      sum += cps->cp - first;
      cpct ++;
    } else {
      bd_tmp = (BD_t *) malloc(sizeof(BD_t));
      bd_tmp->bd = sum/cpct + first;
      bd_tmp->support = cpct;
      bd_tmp->next = bds;
      bds = bd_tmp;
      first = cps->cp;
      sum = 0;
      cpct = 1;
    }
    if (cps->next == NULL) {
      bd_tmp = (BD_t *) malloc(sizeof(BD_t));
      bd_tmp->bd = sum/cpct + first;
      bd_tmp->support = cpct;
      bd_tmp->next = bds;
      bds = bd_tmp;
    }
    last = cps->cp;
    cps = cps->next;
  }
  return bds;
}
#endif



BD_t *CP_cluster(CP_t *cps) {
  int32_t first = cps->cp, last = cps->cp, sum = 0;
  CP_t *begin = cps;
  int cpct = 0;
  BD_t *bds = NULL, *bd_tmp;

  while (cps != NULL) {
    if (cps->cp - first <= 20 && cps->cp - last <= 10) {
      sum += cps->cp - first;
      cpct ++;
      last = cps->cp;
      cps = cps->next;
    } else {
      bd_tmp = (BD_t *) malloc(sizeof(BD_t));
      bd_tmp->bd = sum/cpct + first;
      bd_tmp->support = cpct;
      bd_tmp->next = bds;
      bds = bd_tmp;
      if (cps->cp - last <= 10)	begin = cps;
      else begin = begin->next;
      if (begin) {
	first = begin->cp;
	last = first;
      }
      sum = 0;
      cpct = 0;
      cps = begin;
    }
  }

  bd_tmp = (BD_t *) malloc(sizeof(BD_t));
  bd_tmp->bd = sum/cpct + first;
  bd_tmp->support = cpct;
  bd_tmp->next = bds;
  bds = bd_tmp;

  return bds;
}



void CP_sort(CP_t **CP_ptr) {
  int i=0, j;
  CP_t *cur, **CPs;

  cur = *CP_ptr;
  while (cur != NULL) {
    i ++;
    cur = cur->next;
  }
  CPs = (CP_t **) malloc(i*sizeof(CP_t *));
  i = 0;
  cur = *CP_ptr;
  while (cur != NULL) {
    *(CPs+i) = cur;
    i ++;
    cur = cur->next;
  }
  qsort(CPs, i, sizeof(CP_t *), CP_cmp);
  i --;
  for (j=0; j<i; j++) {
    (*(CPs+j))->next = *(CPs+j+1);
  }
  (*(CPs+i))->next = NULL;
  *CP_ptr = *CPs;

  free(CPs);
}




void BD_sort(BD_t **BD_ptr) {
  int i=0, j;
  BD_t *cur, **BDs;

  cur = *BD_ptr;
  while (cur != NULL) {
    i ++;
    cur = cur->next;
  }

  BDs = (BD_t **) malloc(i*sizeof(BD_t*));
  i = 0;
  cur = *BD_ptr;;
  while (cur != NULL) {
    *(BDs+i) = cur;
    i ++;
    cur = cur->next;
  }

  qsort(BDs, i, sizeof(BD_t *), BD_cmp);

  i --;
  for (j=0; j<i; j++) {
    (*(BDs+j))->next = *(BDs+j+1);
  }
  (*(BDs+i))->next = NULL;;
  *BD_ptr = *BDs;

  free(BDs);
}




int span(ELEMENT_t *ele, int32_t cut) {
  /* span requires PCP sorted according to cp_cmp and images sorted according to frag_cmp */
  int left=0, rite=0;
  IMG_DATA_t *id;

  id = ele->to_img_data;
  while (id && id->to_image->frag.lb <= cut-10) {
    left ++;
    if (id->to_image->frag.rb <= cut+10) rite ++;
    id = id->next;
  }
  return (left-rite)*FUDGE;
}





void TBD_merge(ELEMENT_t *ele) {
  BD_t *prev=NULL, *cur, *next;

  BD_sort(&ele->TBD);
  cur = ele->TBD;
  next = cur->next;
  while (next != NULL) {
    if (next->bd - cur->bd <= 10) {
      if (cur->support < next->support) {
	if (prev == NULL) {
	  ele->TBD = next;
	} else {
	  prev->next = next;
	}
	free(cur);
	cur = next;
      } else {
	cur->next = next->next;
	free(next);
      }
      if (cur) next = cur->next;
      else next = NULL;
    } else {
      prev = cur;
      cur = next;
      next = cur->next;
    }
  }
}





/**************************
 * Dissecting the element *
 **************************/



void dissect(ELE_INFO_t *ele_info) {
  IMG_DATA_t *cur_img_data, *img_data_tmp, *next;
  MSP_t *msp_tmp, *msp_ori;
  IMAGE_t *img_partner, *target_img, *target_partner;
  ELEMENT_t *ele_partner, *cur_ele = ele_info->ele;
  BD_t *tbd_tmp, *new_tbd;
  short dissected;

    cur_img_data = cur_ele->to_img_data;
    while (cur_img_data != NULL) {
      next = cur_img_data->next;
      dissected = 0;
      img_partner = partner(cur_img_data->to_image);
      ele_partner = img_partner->ele_info->ele;
      if (full_length(img_partner, CUTOFF2)){
        ele_partner->flimg_no --;
      }
      tbd_tmp = cur_ele->TBD;
      while (tbd_tmp != NULL) {
	if (tbd_tmp->bd > cur_img_data->to_image->frag.lb || !tbd_tmp->next) {
	  if (tbd_tmp->bd > cur_img_data->to_image->frag.lb && tbd_tmp->bd < cur_img_data->to_image->frag.rb) {
	    if (tbd_tmp->bd - cur_img_data->to_image->frag.lb <= TOO_SHORT) {
	      cur_img_data->to_image->to_msp->score = (int32_t) (cur_img_data->to_image->frag.rb-tbd_tmp->bd+1.)/(cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.)*cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      else img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd;
	    } else if (tbd_tmp->bd - cur_img_data->to_image->frag.rb >= -TOO_SHORT) {
	      cur_img_data->to_image->to_msp->score = (int32_t) (tbd_tmp->bd-cur_img_data->to_image->frag.lb+1.)/(cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.)*cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.rb -= cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      else img_partner->frag.lb += cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      cur_img_data->to_image->frag.rb = tbd_tmp->bd;
	    } else {
	      dissected = 1;
	      /* create a new MSP to hold the left products of the dissection */
	      msp_tmp = add_msp(cur_img_data->to_image->to_msp);
	      /* upgrade the content of the new MSP */
	      if (cur_img_data->to_image == &cur_img_data->to_image->to_msp->query) { 
		target_img = &msp_tmp->query;
		target_partner = &msp_tmp->sbjct;
	      } else {
		target_img = &msp_tmp->sbjct;
		target_partner = &msp_tmp->query;
	      }
	      msp_tmp->score = (int32_t) (tbd_tmp->bd-target_img->frag.lb+1.)/(target_img->frag.rb-target_img->frag.lb+1.)*msp_tmp->score;
	      if (msp_tmp->direction == 1) {
                // RMH: When tracking down the divide by zero problem I ended up signaling this 
                //      spot in the code as the likely place where things went initially wrong. 
                //      Explore further.
		target_partner->frag.rb -= target_img->frag.rb - tbd_tmp->bd;
	      }
	      else {
		target_partner->frag.lb += target_img->frag.rb - tbd_tmp->bd;
	      }
	      target_img->frag.rb = tbd_tmp->bd;
	      fprint_msp(new_msps, msp_tmp);
	      /* keep or ignore the new msp */
	      register_image(target_img, cur_ele);

	      /* create the right product */
	      /* we cheat here.  instead of generating a new MSP, we change the content of the original one, then output it into new_msp when finished.  in the meantime, the index of the MSP is unchanged, which points to the oringinal MSP */
	      cur_img_data->to_image->to_msp->score -= msp_tmp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      else img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd +1;
	    }
	  }
	  if (tbd_tmp->bd >= cur_img_data->to_image->frag.rb || !tbd_tmp->next) { /* end of this image */
	    if (dissected) {
		fprint_msp(new_msps, cur_img_data->to_image->to_msp);
	    }
	    if (too_short(&cur_img_data->to_image->frag) || too_short(&img_partner->frag)) {
	      if (next && next->to_image->to_msp == cur_img_data->to_image->to_msp) next = next->next;
	      dump_image(cur_img_data->to_image);
	    }
	    break;
	  }
	}
	tbd_tmp = tbd_tmp->next;
      }
      cur_img_data = next;
    }
    /* the reason for not generating new elements according to the TBDs and then dissect and partition the images is that there's no gaining in time, it's virtually the same thing */
}




int too_short(FRAG_t *f) {
  if (f->rb - f->lb <= TOO_SHORT) return 1;
  return 0;
}


MSP_t *add_msp(MSP_t *m) {
    MSP_t *msp_tmp;
    /*MSP_DATA_t *hanger;*/

    msp_tmp = MSP_malloc();
    /*hanger = msp_tmp->hanger;*/
    *msp_tmp = *m;
    /*msp_tmp->hanger = hanger;*/
    msp_index ++;
    msp_tmp->query.to_msp = msp_tmp;
    msp_tmp->query.index = 2*msp_index;
    msp_tmp->sbjct.to_msp = msp_tmp;
    msp_tmp->sbjct.index = 2*msp_index+1;
    msp_tmp->stat = 's';
    /*msp_in_mem ++;*/
    /*msp_ct ++;*/

    return msp_tmp;
}




void register_image(IMAGE_t *i, ELEMENT_t *ele) {
  IMAGE_t *ip = partner(i);
  int ct;

  if (too_short(&i->frag)) {
    /* rec_image(&i->ele_info->ele->ignored, i);
    rec_image(&ip->ele_info->ele->ignored, ip); */
    /* the msp here havn't been added to the corresponding elements, so we can simply free it */
    MSP_free(i->to_msp);
  } else {
      /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
      if (ct != i->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  image tree changed before inserting a new image: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
      }*/
      insert_image(&i->ele_info->ele->to_img_tree, i);
      i->ele_info->ele->img_no ++;
      /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
      if (ct != i->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble inserting a new image: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
      }*/
      if (i->ele_info->ele->to_img_data) put_image(&i->ele_info->ele->to_img_data, i);
      /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
      if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  tree changed before inserting a new image: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
      }*/
      insert_image(&ip->ele_info->ele->to_img_tree, ip);
      ip->ele_info->ele->img_no ++;
      /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);      
      if (ct != ip->ele_info->ele->img_no) {      
	err_no ++;
        fprintf(log_file, "error:  trouble inserting a new image: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct); 
	fflush(log_file);
        exit(2);      
      }*/      
      if (ip->ele_info->ele->to_img_data) put_image(&ip->ele_info->ele->to_img_data, ip);
  }
}





void put_image(IMG_DATA_t **idp, IMAGE_t *i) {
  IMG_DATA_t *img_data_tmp;

  img_data_tmp = (IMG_DATA_t *)malloc(sizeof(IMG_DATA_t));
  img_data_tmp->to_image = i;
  img_data_tmp->next = *idp;
  *idp = img_data_tmp;
}




void dump_image(IMAGE_t *i) {
    IMAGE_t *ip = partner(i);
    int ct;

    if (i->ele_info->ele->to_img_data) remove_image(i);
    if (ip->ele_info->ele->to_img_data) remove_image(ip);

    /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
    if (ct != i->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  image tree changed before deleting a node: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
    }*/
    delete_image(&i->ele_info->ele->to_img_tree, i);
    i->ele_info->ele->img_no --;
    /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
    if (ct != i->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble deleting an image node: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/

    /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
    if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  image tree changed before deleting a node: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/
    delete_image(&ip->ele_info->ele->to_img_tree, ip);
    ip->ele_info->ele->img_no --;
    /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
    if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble deleting an image node: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/
    
    MSP_free(i->to_msp);
}




void remove_image(IMAGE_t *i) {
    IMG_DATA_t *prev_img_data=NULL, *cur_img_data;

    cur_img_data = i->ele_info->ele->to_img_data;
    while (cur_img_data != NULL) {
      if (cur_img_data->to_image == i) {
	if (prev_img_data != NULL) {
	  prev_img_data->next = cur_img_data->next;
	} else {
	  i->ele_info->ele->to_img_data = cur_img_data->next;
	}
	free(cur_img_data);
	break;
      }
      prev_img_data = cur_img_data;
      cur_img_data = cur_img_data->next;
    }
}




/*****************************
 *updating dissected element *
 *****************************/





void combo_update(ELE_INFO_t *ele_info) {
  if (ele_info->ele->img_no < 0) {
    err_no ++;
    fprintf(log_file, "error:  combo ele %d has %d images\n", ele_info->index, ele_info->ele->img_no);
    fflush(log_file);
    exit(2);
  }
    ele_info->stat = 'X';
    if (ele_info->ele->edges) combo_edge_update(ele_info, &ele_info->ele->edges);
    if (ele_info->ele->edge_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d edge_node left\n", ele_info->index, ele_info->ele->edge_no);
	fflush(log_file);
	exit(4);
    }
    ele_info->ele->flimg_no = 0;
    if (ele_info->ele->PCP) CP_free(&ele_info->ele->PCP);
    if (ele_info->ele->redef) {
      if (ele_info->ele->to_img_data) {
	err_no ++;
	fprintf(log_file, "error re-defining ele %d, still images left\n", ele_info->index);
	fflush(log_file);
	exit(5);
      }
      /* all images and msps are relocated to offsprings, so img_tree_free() is enough */
      if (ele_info->ele->to_img_tree) img_tree_free(&ele_info->ele->to_img_tree, ele_info);
      if (ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d img_node left\n", ele_info->index, ele_info->ele->img_no);
	fflush(log_file);
	ele_info->ele->img_no = 0;
	exit(2);
      }
      combo_output(ele_info);
    } else {
      if (ele_info->ele->img_no || ele_info->ele->to_img_data || ele_info->ele->to_img_tree) {
	err_no ++;
	fprintf(log_file, "error:  images not cleaned in obs ele %d\n", ele_info->index);
	fflush(log_file);
	exit(2);
      }
      else obs_output(ele_info);
    }
}





void combo_edge_update(ELE_INFO_t *ele_info, EDGE_TREE_t **edge_node_p) {
  ELE_INFO_t *epi;
  EDGE_TREE_t *edge_node = *edge_node_p;
  ELEMENT_t *ele = ele_info->ele;

  if (edge_node->l) combo_edge_update(ele_info, &(*edge_node_p)->l);

  if (edge_node->r) combo_edge_update(ele_info, &(*edge_node_p)->r);

  epi = linked_ele(ele_info, edge_node->to_edge);
  if (epi->stat == 'v') epi->stat = 'w';
  if (epi->ele->PCP) CP_clean(&epi->ele->PCP, ele_info);
  delete_edge(&epi->ele->edges, edge_node->to_edge);
  epi->ele->edge_no --;

  EDGE_free(edge_node->to_edge);

  ele_info->ele->edge_no --;
  free(edge_node);
  *edge_node_p = NULL;
}





void CP_clean(CP_t **cps_p, ELE_INFO_t *cont) {
  /*int res=0;*/
  CP_t *cp_cur, *cp_prev;

  cp_cur = *cps_p;
  cp_prev = NULL;
  while (cp_cur != NULL) {
    if (cp_cur->contributor == cont) {
      if (cp_prev == NULL) {
	*cps_p = cp_cur->next;
	free(cp_cur);
	cp_cur = *cps_p;
      } else {
	cp_prev->next = cp_cur->next;
	free(cp_cur);
	cp_cur = cp_prev->next;
      }
      /*res ++;*/
    } else {
      cp_prev = cp_cur;
      cp_cur = cp_cur->next;
    }
  }
  /*return res;*/
}





/****************************************************
 ****************************************************
 **  functions for finding PCPs and set up edges   **
 ****************************************************
 ****************************************************/



void edges_and_cps(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) { /* the core function in this part, calls everything below */
    int eff_img_ct, i, prim;
    int j, riteplace;
    short ritetime;
    MSP_t *prim_p;

    IMG_DATA_t *cur_img_data;
    ELEMENT_t *ele_partner, *cur_ele = ele_info->ele;
    ELE_INFO_t *epi;
    IMAGE_t *img_partner, *cur_img;

    IMAGE_t *token_image;
    int token_mark=0;
    IMG_NODE_t *consis_rt, *cur_consis_nd;
    ELE_DATA_t *ele_data_tmp;
    IMG_DATA_t *scovs;

    int32_t max_score=0;
    short dir;

    /* sort unprocessed images according to their partner element, and then their left bounds */
    /* to allocate proper amount of memory */
    /* when update is 1, it means the element is an offspring of a combo, in which case all the images need to be processed; when update is 0, those images whose partner is 'v' or 'w' should be omitted */
    /* no self-images are accepted to generate pcps or edges, updating or not */
    if (!cur_ele->to_img_data) listify(cur_ele->to_img_tree, &cur_ele->to_img_data);
    cur_img_data = cur_ele->to_img_data;
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info;
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') eff_img_ct ++;
      }
      cur_img_data = cur_img_data->next;
    }

  if (eff_img_ct) {
    token_image = (IMAGE_t *) malloc(sizeof(IMAGE_t));
    token_image->frag.lb = 0;
    token_image->frag.rb = 0;
    token_image->to_msp = NULL;
    token_image->ele_info = NULL;
    consis_rt = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
    consis_rt->to_image = token_image;
    consis_rt->children = NULL;
    consis_rt->sib = NULL;

    if (eff_img_ct > MAX_IMG) {
      img_ptr = (IMAGE_t **) malloc(eff_img_ct*sizeof(IMAGE_t *));
    }

    cur_img_data = cur_ele->to_img_data;
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info; 
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') {
  	  *(img_ptr+eff_img_ct) = cur_img_data->to_image;
	  eff_img_ct ++;
        }
      }
      cur_img_data = cur_img_data->next;
    }
    qsort(img_ptr, eff_img_ct, sizeof(IMAGE_t *), partner_cmp);

    /* recognizing full-length images, and put partial and secondary images into a consistency tree, in which images connected are consistent with each other (look at the consis() function of definition of consistency) */
    ritetime = 1;
    for (i=0; i<eff_img_ct; i++) {
      cur_img = *(img_ptr+i);
      img_partner = partner(cur_img);
      if (ritetime) { /* new ele_partner begins */
	epi = img_partner->ele_info;
	ritetime = 0;
	if (!epi->ele) ele_read_in(epi, 1);
	ele_partner = epi->ele;
	riteplace = i;
	prim = 0;
	prim_p = NULL;
	max_score = 0;
      }
      if (img_partner->ele_info->index == epi->index) { /* still the same partner element */
	/* full length */
	if (full_length(cur_img, CUTOFF2)) {
	  /*cur_img->to_msp->stat = 'p';*/
	  prim = 1;
	  cur_ele->flimg_no ++;
	}
	if (full_length(img_partner, CUTOFF2)) {
	  /*cur_img->to_msp->stat = 'p';*/
	  prim = 1;
	  ele_partner->flimg_no ++;
	}
	/*if (cur_img->to_msp->stat == 'p') {*/
	if (prim == 1) {
	  prim = 0;
	  if (cur_img->to_msp->iden > max_score) {
	    max_score = cur_img->to_msp->iden;
	    dir = cur_img->to_msp->direction;
	    prim_p = cur_img->to_msp;
	  }
	}
	if (!prim_p) { /* non-fl, and no fl found yet */
	  consis_tree_build(consis_rt, cur_img, 1);
	}
      }
      if (img_partner->ele_info->index != epi->index || i == eff_img_ct-1) { /* reached the end of the current ele_partner */
	/* start a new ele_partner */
	ritetime = 1;
	if (img_partner->ele_info->index != epi->index) i --;
	/* finish the current ele_partner */
	if (prim_p) {
	  prim_p->stat = 'p';
	  prim = 1;
	} else { /* if no fl, find partial primary images */
	  /*prim = find_prim(consis_rt, CUTOFF2, 0, 0, 0, &max_score, &dir);*/
	  prim = find_prim(consis_rt->children, CUTOFF2, ele_info->ele->frag.lb, -1, 0, 0, 0, 0, 0, &token_mark, &max_score, &dir);
	}
	/* build edge */
	if (prim) {
	  if (ele_info->index != epi->index) add_edge(ele_info, epi, 'p', max_score, dir);
	  else {
	    err_no ++;
	    fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
	    fflush(log_file);
	  }
	  for (j=riteplace; j<=i; j++) {
	    if ((*(img_ptr+j))->to_msp->stat == 'p') {
	      cur_img = *(img_ptr+j);
	      img_partner = partner(cur_img);
	      add_CP(&cur_ele->PCP, cur_img->frag.lb, epi);
	      add_CP(&cur_ele->PCP, cur_img->frag.rb, epi);
	      add_CP(&ele_partner->PCP, img_partner->frag.lb, ele_info);
	      add_CP(&ele_partner->PCP, img_partner->frag.rb, ele_info);
	    }
	  }
	} else { /* no primary images, full-length or partial */
	  if (ele_info->index != epi->index) add_edge(ele_info, epi, 's', 0, 0);
          else {
            err_no ++;
            fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
	    fflush(log_file);
          }
	}
	if (consis_rt->children != NULL) {
	  consis_tree_free(consis_rt->children);
	  consis_rt->children = NULL;
	}
      }
    }

    if (eff_img_ct > MAX_IMG) {
      free(img_ptr);
    }
    free(consis_rt);
    free(token_image);
  }

  if (ele_info->ele->edges) adjust_edge_tree(ele_info); 

  cur_ele->update = 0;
  ele_info->stat = 't';
}



/* NOTE: Due to the difference in computational precesion between 32bit
 *       and 64 bit processors the following function may return differing
 *       results.
 */
int full_length(IMAGE_t *i, float cutoff) {
  if (!i->ele_info->ele) {
    err_no ++;
    fprintf(log_file, "error:  element %d not in memory\n", i->ele_info->index);
    fflush(log_file);
    exit(3);
    ele_read_in(i->ele_info, 1);
  }
  if (i->frag.lb-i->ele_info->ele->frag.lb < 10 && i->frag.rb-i->ele_info->ele->frag.rb > -10 && ((float) i->frag.rb-i->frag.lb)/(i->ele_info->ele->frag.rb-i->ele_info->ele->frag.lb) > cutoff) {
    return 1;
  }
  return 0;
}



void add_CP(CP_t **CP_ptr, int32_t cp, ELE_INFO_t *cont) {
  CP_t *new = (CP_t *) malloc(sizeof(CP_t));
  new->cp = cp;
  new->contributor = cont;
  new->next = *CP_ptr;
  *CP_ptr = new;
}



void add_edge(ELE_INFO_t *ele1_info, ELE_INFO_t *ele2_info, char type, int32_t score, short dir) {
  EDGE_t *new = EDGE_malloc();
  int ct;

  edge_index ++;
  new->index = edge_index;
  new->ele1_info = ele1_info;
  new->ele2_info = ele2_info;
  new->type = type;
  new->score = score;
  new->direction = dir;
  /*ct = count_edge_nodes(ele1_info->ele->edges);
  if (ct != ele1_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble before inserting edge\n");
    fflush(log_file);
    exit(4);
  }*/
  ele1_info->ele->edge_no ++;
  insert_edge(&ele1_info->ele->edges, new);
  /*ct = count_edge_nodes(ele1_info->ele->edges);
  if (ct != ele1_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble after inserting edge\n");
    fflush(log_file);
    exit(4);
  }*/
  if (ele1_info->index != ele2_info->index) {
    /*ct = count_edge_nodes(ele2_info->ele->edges);
    if (ct != ele2_info->ele->edge_no) {
      err_no ++;
      fprintf(log_file, "error:  trouble before inserting edge\n");
      fflush(log_file);
      exit(4);
    }*/
    ele2_info->ele->edge_no ++;
    insert_edge(&ele2_info->ele->edges, new);
    /*ct = count_edge_nodes(ele2_info->ele->edges);    
    if (ct != ele2_info->ele->edge_no) {    
      err_no ++;    
      fprintf(log_file, "error:  trouble after inserting edge\n");    
      fflush(log_file);
      exit(4);
    }*/  
  }
}




void adjust_edge_tree(ELE_INFO_t *ele_info) {
  EDGE_t **edge_array;
  int ct, ct1;

  /*ct = count_edge_nodes(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble before adjusting the edge tree: ele %d %d %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  }
  ct = count_total_edges(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  illegitimate edge in ele %d\n", ele_info->index);
    fflush(log_file);
    exit(5);
  }*/
  ct = ele_info->ele->edge_no;
  edge_array = (EDGE_t **) malloc(ct*sizeof(EDGE_t *));
  ct1 = charge_edge_array(edge_array, ele_info->ele->edges, 0);
  if (ct1 != ct) {
    err_no ++;
    fprintf(log_file, "error:  trouble charging the edge array in ele %d: %d charged, %d expected\n", ele_info->index, ct1, ct);
    fflush(log_file);
    exit(4);
  }
  edge_tree_free(&ele_info->ele->edges);
  build_edge_tree(&ele_info->ele->edges, edge_array, 0, ct-1);
  /*ct = count_edge_nodes(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble adjusting the edge tree: ele %d %d %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  }*/
  free(edge_array);
}




int charge_edge_array(EDGE_t **edge_array, EDGE_TREE_t *rt, int pos) {
  if (rt->l) pos = charge_edge_array(edge_array, rt->l, pos);

  *(edge_array+pos) = rt->to_edge;
  pos ++;

  if (rt->r) pos = charge_edge_array(edge_array, rt->r, pos);

  return pos;
}




/******************************************
 * functions for building the consis_tree *
 ******************************************/


int consis_tree_build(IMG_NODE_t *rt, IMAGE_t *im, int prequal) {
  int sum=0;
  IMG_NODE_t *nex_rt, *node;
  if (prequal || consis(rt->to_image, im, CUTOFF2)) {
    nex_rt = rt->children;
    while (nex_rt != NULL) {
      sum += consis_tree_build(nex_rt, im, 0);
      nex_rt = nex_rt->sib;
    }
    if (!sum) {
      node = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
      node->recorded = 0;
      node->to_image = im;
      node->sib = NULL;
      node->children = NULL;
      *node_entry(&rt->children) = node;
      sum = 1;
    }
  }
  return sum;
}


int consis(IMAGE_t *i1, IMAGE_t *i2, float cutoff) {
  int res = 0;
  IMAGE_t *ip1 = partner(i1), *ip2 = partner(i2);
  if (i1->ele_info->index == i2->ele_info->index && ip1->ele_info->index == ip2->ele_info->index && i1->to_msp->direction == i2->to_msp->direction) {
    if (i1->to_msp->direction == 1) {
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.lb - ip2->frag.lb) > 0) {
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) && !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
	  res = 1;
	}
      }
    } else {
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.rb - ip2->frag.rb) < 0) {
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) && !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
	  res = 1;
	}
      }
    }
  }
  return res;
}



IMG_NODE_t ** node_entry(IMG_NODE_t **node_pp) {
  if (*node_pp != NULL) {
    return node_entry(&(*node_pp)->sib);
  }
  return node_pp;
}




void consis_tree_free(IMG_NODE_t *rt) {
  if (rt->sib != NULL) {
    consis_tree_free(rt->sib);
  } 
  if (rt->children != NULL) {
    consis_tree_free(rt->children);
  }
  free(rt);
}




/*****************************************
 * functions for parsing the consis_tree *
 *****************************************/

/*
  I. function for finding primary images from the consis_tree
*/

#if 0
int find_prim(IMG_NODE_t *nd, float cutoff, int32_t score, int32_t hist, short which, int32_t *sc, short *d) {
  IMG_NODE_t *nex_node;
  int sum = 0;
  short level, further = 0;
  IMAGE_t *img_partner;

  if (!hist) level = 1;
  hist += nd->to_image->frag.rb - nd->to_image->frag.lb;
  
  if (level && hist) { /* if first image in a path isn't EOE, forget the path */
    if (nd->to_image->frag.lb - nd->to_image->ele_info->ele->frag.lb < MARGIN) {
      which = 1;
    }
    img_partner = partner(nd->to_image);
    if (nd->to_image->to_msp->direction == 1) {
      if (img_partner->frag.lb - img_partner->ele_info->ele->frag.lb < MARGIN) {
	which += 2;
      }
    } else {
      if (img_partner->frag.rb - img_partner->ele_info->ele->frag.rb > -MARGIN) {
	which += 2;
      }
    }
    if (!which) return 0;
  }

  if (nd->children != NULL) {
    nex_node = nd->children;
    while (nex_node != NULL) {
      sum += find_prim(nex_node, cutoff, score, hist, which, sc, d);
      nex_node = nex_node->sib;
    }
    if (sum && hist) {
      nd->to_image->to_msp->stat = 'p'; /* p stands for primary, NOT partial */
    }
  } else { /* the last image goes to EOE on the same element as the first image in the path */
    if ((which == 1 || which == 3) && nd->to_image->frag.rb - nd->to_image->ele_info->ele->frag.rb > -MARGIN  && (float) hist/(nd->to_image->ele_info->ele->frag.rb - nd->to_image->ele_info->ele->frag.lb) > cutoff) {
      further = 1;
    } else {
      img_partner = partner(nd->to_image);
      if ((which == 2 || which == 3) && (float) hist/(img_partner->ele_info->ele->frag.rb - img_partner->ele_info->ele->frag.lb) > cutoff) {
	if (nd->to_image->to_msp->direction == 1) {
	  if (img_partner->frag.rb - img_partner->ele_info->ele->frag.rb > -MARGIN) further = 1;
	} else {
	  if (img_partner->frag.lb - img_partner->ele_info->ele->frag.lb < MARGIN) further = 1;
	}
      }
    }
    if (further) {
      nd->to_image->to_msp->stat = 'p';
      score += nd->to_image->to_msp->score;
      if (score > *sc) {
	*sc = score;
	*d = nd->to_image->to_msp->direction;
      }
      sum = 1;
    } else {sum = 0;}
  }

  return sum;
}
#endif



int find_prim(IMG_NODE_t *nd, float cutoff, int32_t end1, int32_t end2, int32_t efl1, int32_t efl2, int32_t al1, int32_t al2, int32_t score, int *pmarkp, int32_t *sc, short *d) {
  int sum = 0, mark=0;
  int32_t skip1, skip2, len1, len2;
  IMAGE_t *ipt;

  if (nd->sib) sum += find_prim(nd->sib, cutoff, end1, end2, efl1, efl2, al1, al2, score, pmarkp, sc, d);

  ipt = partner(nd->to_image);
  if (end2 < 0) { /*first alignment in the group*/
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->ele_info->ele->frag.lb;
    else end2 = ipt->ele_info->ele->frag.rb;
  }

  skip1 = nd->to_image->frag.lb - end1;
  if (nd->to_image->to_msp->direction == 1) skip2 = ipt->frag.lb - end2;
  else skip2 = end2 - ipt->frag.rb;
  if (skip1>10 && skip2>10) {
    efl1 += skip1;
    efl2 += skip2;
  }
  len1 = nd->to_image->frag.rb - nd->to_image->frag.lb;
  len2 = ipt->frag.rb - ipt->frag.lb;
  efl1 += len1;
  efl2 += len2;
  al1 += len1;
  al2 += len2;
  score += ((int32_t) nd->to_image->to_msp->iden)*(len1+len2)/2;

  if (nd->children) {
    end1 = nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->frag.rb;
    else end2 = ipt->frag.lb;
    sum += find_prim(nd->children, cutoff, end1, end2, efl1, efl2, al1, al2, score, &mark, sc, d);
  } else { /*last alignment in group*/
    skip1 = nd->to_image->ele_info->ele->frag.rb - nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1) skip2 = ipt->ele_info->ele->frag.rb - ipt->frag.rb;
    else skip2 = ipt->frag.lb - ipt->ele_info->ele->frag.lb;
    if (skip1>10 && skip2>10) {
      efl1 += skip1;
      efl2 += skip2;
    }
    /*if (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) {*/
    if ( (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) && (efl1-al1 < 30 || efl2-al2 < 30) ) {
      sum = 1;
      mark = 1;
      if ( (al1+al2) == 0 )
      {
         // RMH: A divide by zero error has occured a few times at this point. 
         //      I looked extensively through the code and cannot determine what
         //      is causing this -- although I have a reliable test case to 
         //      explore this further.  For now this is the only workaround.
         //      If al1 & al2 are zero then simply set them to 1.  
         printf("eleredef warning: Divide by zero averted -- setting al1 to 1.\n");
         al1 = 1;
      }
      score = score*2/(al1+al2);
      if (score > *sc) {
	*sc = score;
	*d = nd->to_image->to_msp->direction;
      }
    }
  }
  if (mark) {
      nd->to_image->to_msp->stat  = 'p';
      *pmarkp = 1;
  }
  return sum;
}






/************************
 ************************
 ***                  ***
 ***      OUTPUT      ***
 ***                  ***
 ************************
 ************************/





void combo_output(ELE_INFO_t *ele_info) {
  char *ele_name;
  FILE *ele_fp;
  ELE_DATA_t *cur_ele_data;

  fprintf(combo, "%d %s %d %d \n", ele_info->index, ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(combo);

#if 0
  /* no ele_read_in() needed, 'cuz where this function is called, both the combo and its offsprings are in the memory */
  ele_name = (char *) malloc(30*sizeof(char));
  sprintf(ele_name, "tmp2/clan/combos/e%d", ele_info->index);
  ele_fp = fopen(ele_name, "w");

  fprintf(ele_fp, "ele %d\n", ele_info->index);
  fprintf(ele_fp, "offsprings \n");
  cur_ele_data = ele_info->ele->redef;
  while (cur_ele_data != NULL) {
    if (cur_ele_data->ele_info->ele->direction == 1) fprintf(ele_fp, "%6d %10s %8d %8d\n", cur_ele_data->ele_info->index, cur_ele_data->ele_info->ele->frag.seq_name, cur_ele_data->ele_info->ele->frag.lb, cur_ele_data->ele_info->ele->frag.rb);
    else fprintf(ele_fp, "%6d %10s %8d %8d\n", cur_ele_data->ele_info->index, cur_ele_data->ele_info->ele->frag.seq_name, cur_ele_data->ele_info->ele->frag.rb, cur_ele_data->ele_info->ele->frag.lb);
    cur_ele_data = cur_ele_data->next;
  }

  fclose(ele_fp);
  free(ele_name);
  /*  remove_ele(ele_info);*/
#endif
}





void obs_output(ELE_INFO_t *ele_info) {
  char *ele_name;
  FILE *ele_fp;

  fprintf(obs, "%d %s %d %d \n", ele_info->index, ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(obs);

#if 0
  ele_name = (char *) malloc(30*sizeof(char));
  sprintf(ele_name, "tmp2/clan/obsolete/e%d", ele_info->index);
  ele_fp = fopen(ele_name, "w");
  fprint_ele_obs(ele_fp, ele_info);
  fclose(ele_fp);
  free(ele_name);
  /*  remove_ele(ele_info);*/
#endif
}




void fprint_ele_obs(FILE *fp, ELE_INFO_t *ele_info) {
  int i;
  BD_t *cur_bd;
  ELEMENT_t *ele=ele_info->ele;

  fprintf(fp, "ele %d\n", ele_info->index);
  fprintf(fp, "%s %d %d \n", ele->frag.seq_name, ele->frag.lb,  ele->frag.rb);
  if (ele->TBD != NULL) {
    fprintf(fp, "cutting points\n");
    i = 0;
    cur_bd = ele->TBD;
    while (cur_bd != NULL) {
      fprintf(fp, "%d ", cur_bd->bd);
      i ++;
      if (i%10 == 0) fprintf(fp, "\n");
      cur_bd = cur_bd->next;
    }
  }
}






/**************************************
 **************************************
 *                                    *
 *  comparison fucntions for qsort()  *
 *                                    *
 **************************************
 **************************************/






int frag_cmp(const void *i1, const void *i2) {
  int res = (*(IMG_DATA_t **)i1)->to_image->frag.seq_name - (*(IMG_DATA_t **)i2)->to_image->frag.seq_name; /*strncmp((*(IMG_DATA_t **)i1)->to_image->frag.seq_name, (*(IMG_DATA_t **)i2)->to_image->frag.seq_name, NAME_LEN)*/
  if (res == 0) {
    res = (*(IMG_DATA_t **)i1)->to_image->frag.lb - (*(IMG_DATA_t **)i2)->to_image->frag.lb;
    if (res == 0) {
      res = (*(IMG_DATA_t **)i1)->to_image->frag.rb - (*(IMG_DATA_t **)i2)->to_image->frag.rb;
    }
  }
  return res;
}


int partner_cmp(const void *i1, const void *i2) {
  int res = partner(*((IMAGE_t **)i1))->ele_info->index - partner(*((IMAGE_t **)i2))->ele_info->index;
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->to_msp->direction - (*((IMAGE_t **)i2))->to_msp->direction;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.lb - (*((IMAGE_t **)i2))->frag.lb;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.rb - (*((IMAGE_t **)i2))->frag.rb;
  }
  return res;
}


int CP_cmp(const void *cp1, const void *cp2) {
  return (*((CP_t **) cp1))->cp - (*((CP_t **) cp2))->cp;
}


int BD_cmp(const void *bd1, const void *bd2) {
  return (*((BD_t **) bd1))->bd - (*((BD_t **) bd2))->bd;
}


int fam_cmp(const void *fd1, const void *fd2) {
  return (*((FAM_DATA_t **) fd1))->to_family->index - (*((FAM_DATA_t **) fd2))->to_family->index;
}
