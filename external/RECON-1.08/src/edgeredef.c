#include "ele.h"

#define CUTOFF3 0.7




void report_redef_stat();

void general_edge_redef(ELE_INFO_t *);
void build_local_edge_net(ELE_INFO_t *, int, ELE_DATA_t **);
void local_edge_net_walk(ELE_INFO_t *, EDGE_TREE_t *, int, ELE_DATA_t **);
void cruise_local_edge_net(ELE_DATA_t *);
void dissolve_local_edge_net(ELE_DATA_t **);
void stage_exit(int );
void edge_filt(ELE_INFO_t *);
EDGE_TREE_t *find_PPS(ELE_INFO_t *, EDGE_TREE_t *);
EDGE_TREE_t * edge_update(ELE_INFO_t *, EDGE_TREE_t *, EDGE_TREE_t *);

void edge_repair(ELE_INFO_t *);
EDGE_TREE_t *best_link(ELE_INFO_t *, EDGE_TREE_t *);




int main (int argc, char *argv[]) {
  int i, rounds=0, ei, start;
  short to_march;
  char line[35], stat;
  FILE *seq_list, *ele_no, *redef_stat;

  /* processing command line */
  if (argc == 1) {
    printf("usage: fam_def seq_list start\n where seq_list is the list of sequence names, start is the index of the element to start defining families.  The latter one is optional.\n");
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

  ele_no = fopen("summary/redef_ele_no", "r");
  if(!ele_no) {
    printf("Can not open redef_ele_no.  Exiting\n");
    exit(1);
  }

  redef_stat = fopen("tmp/redef_stat", "r");
  if (!redef_stat) {
    printf("Can not open redef_stat.  Exiting\n");
    exit(1);
  }

  log_file = fopen("tmp2/log", "w");
  if (!log_file) {
    printf("Can not open log_file.  Exiting\n");
    exit(1);
  }

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

  /* update edges, to get rid of PPS */

  for (i=start; i<ele_ct; i++) {
    if ((*(all_ele+i))->stat != 'O' && (*(all_ele+i))->stat != 'X') {
      fprintf(log_file, "evaluating edges of element %d\n", (*(all_ele+i))->index);
      fflush(log_file);
      if ((*(all_ele+i))->stat == 'v') {
	general_edge_redef(*(all_ele+i));
	report_redef_stat();
      } else if((*(all_ele+i))->stat != 'y') {
	err_no ++;
	fprintf(log_file, "ele %d %c not properly defined\n",  (*(all_ele+i))->index, (*(all_ele+i))->stat);
	fflush(log_file);
	exit(5);
      }
    }
  }

  for (i=start; i<ele_ct; i++) {
    if ((*(all_ele+i))->stat != 'O' && (*(all_ele+i))->stat != 'X') {
      fprintf(log_file, "repairing edges of element %d\n", (*(all_ele+i))->index);
      fflush(log_file);
      if ((*(all_ele+i))->stat == 'y') {
        edge_repair(*(all_ele+i));
      } else {
        err_no ++;
        fprintf(log_file, "ele %d %c edges not properly redefined\n",  (*(all_ele+i))->index, (*(all_ele+i))->stat);
        fflush(log_file);
        exit(5);
      }
    }
  }

  fprintf(log_file, "total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  fprintf(log_file, "%d files read, %d msps seen, %d edges seen\n", files_read, msp_ct, edge_ct);
  fprintf(log_file, "%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);


  exit(0);
}




void report_redef_stat() {
  int i;
  FILE *redef_stat;

  redef_stat = fopen("tmp2/redef_stat", "w");
  for (i=0; i<ele_ct; i++) {
    fprintf(redef_stat, "%d %c %d\n", (*(all_ele+i))->index, (*(all_ele+i))->stat, (*(all_ele+i))->file_updated);
  }

  fclose(redef_stat);
}








/***********************************
 ***********************************
 **                               **
 **      UPDATING THE EDGES       **
 **                               **
 ***********************************
 ***********************************/




void general_edge_redef(ELE_INFO_t *ele_info) {
  ELE_DATA_t *local_net=NULL;
  char *command;

  if (!ele_info->ele) ele_read_in(ele_info, 2);

  if (ele_info->ele->edges) {
    clan_ct ++;
    clan_size = 0;
    clan_core_size = 0;
    fprintf(log_file, "new clan: %d for ele %d\n", clan_ct, ele_info->index);
    fflush(log_file);
    build_local_edge_net(ele_info, 1, &local_net); /*notice level here is 0, not 1 */
    fprintf(log_file, "clan size: %d, clan core size: %d\n", clan_size, clan_core_size);
    fflush(log_file);
    cruise_local_edge_net(local_net);
    dissolve_local_edge_net(&local_net);
  } else {
    ele_info->stat = 'y';
    ele_cleanup(&ele_info->ele);
    command = (char *) malloc(80*sizeof(char));
    sprintf(command, "cp tmp/e%d tmp2/.", ele_info->index);
    system(command);
    free(command);
  }
}




/* notice the structural similarity and difference b/t the following two functions and the corresponding ones in redef.c */
/* notice that it's much simpler here */

void build_local_edge_net(ELE_INFO_t *ele_info, int level, ELE_DATA_t **net_p) {
  ELE_INFO_t *epi;
  ELE_DATA_t *que;
  ELEMENT_t *ele;

  clan_size ++;
  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = *net_p;
  *net_p = que;

  if (!ele_info->ele) ele_read_in(ele_info, 1);
  ele = ele_info->ele;

  if (level < DEPTH) {
    if (ele_info->stat != 'y') {
      ele->l_hold = 1;
      clan_core_size ++;
      /*fprintf(log_file, "    %d\n", ele_info->index);
	fflush(log_file);*/
      level ++;
      if (ele->edges) local_edge_net_walk(ele_info, ele->edges, level, net_p);
      else {
	err_no ++;
	fprintf(log_file, "error:  edge tree missing in ele %d\n", ele_info->index);
	fflush(log_file);
      }
    } else ele->l_hold = 2;
  } else ele->l_hold = 2;
}





void local_edge_net_walk(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, int level, ELE_DATA_t **net_p) {
  ELE_INFO_t *epi;
  ELEMENT_t *ele=ele_info->ele;

  if (rt->l) local_edge_net_walk(ele_info, rt->l, level, net_p);

  if (rt->to_edge->type == 'p' || rt->to_edge->type == 'S'  || rt->to_edge->type == 'P') {
    epi = linked_ele(ele_info, rt->to_edge);
    if (!epi->ele) ele_read_in(epi, 2);
    if (!epi->ele->l_hold) build_local_edge_net(epi, level, net_p);
  }

  if (rt->r) local_edge_net_walk(ele_info, rt->r, level, net_p);
}




void cruise_local_edge_net(ELE_DATA_t *local_net) {
  ELE_DATA_t *que;

  que = local_net;
  while (que) {
    /* here, all elements in the local_network are explicitly hung in ques,
       unlike in the local_network for ele_ref, where some are hung in redef
       of other elements */
    /* also, no need to march for a second time, since once edge_filt is done,
       the status will not be affected by the processing of other elements, so
       no need to come back and re-do it */
    if (que->ele_info->ele->l_hold == 1) {
      if (que->ele_info->stat != 'y') {
	if (que->ele_info->ele->edges) edge_filt(que->ele_info);
	else if (que != local_net) {
	  fprintf(log_file, "error: ele %d edge tree missing\n", que->ele_info->index);
	  fflush(log_file);
	  exit(4);
	}
	que->ele_info->stat = 'y';
      }
      else {
	fprintf(log_file, "error:  ele %d conflicting status in local network\n", que->ele_info->index);
	fflush(log_file);
	exit(5);
      }
    }
    que = que->next;
  }
}




void dissolve_local_edge_net(ELE_DATA_t **net_p) {
  ELE_DATA_t *que;
  char *command;
  int i, ele_left=0;

  que = *net_p;
  while (que) {
    que->ele_info->ele->l_hold = 0;
    ele_write_out(que->ele_info, 2);
    ele_cleanup(&que->ele_info->ele);
    que = que->next;
  }
  ele_data_free(net_p);

  if (msp_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d msps total, %d seen, %d left in memory\n", msp_index+1, msp_ct, msp_in_mem);
    fflush(log_file);
    msp_left += msp_in_mem;

  }

  if (edge_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d edges total, %d seen, %d left in memory\n", edge_index+1, edge_ct, edge_in_mem);
    fflush(log_file);
    edge_left += edge_in_mem;
    edge_in_mem = 0;
  }

  for (i=0; i<ele_ct; i++) {
    if ((*(all_ele+i))->ele) {
      ele_left ++;
      err_no ++;
      fprintf(log_file, "error:  element %d not cleaned from memory. stat = %c\n", i+1, (*(all_ele+i))->stat);
      fflush(log_file);
    }
  }
  if (ele_left) {
    err_no ++;
    fprintf(log_file, "error:  %d elements left in memory\n", ele_left);
    fflush(log_file);
  }

  if (err_no) exit(3);

#if 0
  command  = (char *) malloc(50*sizeof(char));
  sprintf(command, "mv -f tmp2/clan/e* tmp2/.");
  if (system(command)) {
    fprintf(log_file, "error moving files from tmp2/clan/\n");
    fflush(log_file);
    stage_exit(6);
  }
  free(command);  
#endif
}




void stage_exit(int code) {
  char *command = (char *) malloc(50*sizeof(char));

  sprintf(command, "mv -f tmp2/redef_stat tmp2/redef_stat_prev");
  if (system(command)) {
    fprintf(log_file, "error in moving tmp2/redef_stat\n");
    fflush(log_file);
    exit(code);
  }

  report_redef_stat();
  exit(code);
}




void edge_filt(ELE_INFO_t *ele_info) {
  EDGE_TREE_t *max;

  max = find_PPS(ele_info, ele_info->ele->edges);

  if (max) {
    max = edge_update(ele_info, ele_info->ele->edges, max);
    max->to_edge->type = 'P';
  }
}



EDGE_TREE_t *find_PPS(ELE_INFO_t *ele_info, EDGE_TREE_t *rt) {
  ELE_INFO_t *epi;
  EDGE_TREE_t *res=NULL;

  if (rt->to_edge->type == 'p' || rt->to_edge->type == 'S') {
    epi = linked_ele(ele_info, rt->to_edge);
    if (((float) ele_info->ele->frag.rb - ele_info->ele->frag.lb)/(epi->ele->frag.rb - epi->ele->frag.lb) < CUTOFF3) res = rt;
  }

  if (!res && rt->l) res = find_PPS(ele_info, rt->l);
  if (!res && rt->r) res = find_PPS(ele_info, rt->r);

  return res;
}




EDGE_TREE_t *edge_update(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, EDGE_TREE_t *ref) {
  ELE_INFO_t *epi;

  if (rt->l) ref = edge_update(ele_info, rt->l, ref);
  if (rt->r) ref = edge_update(ele_info, rt->r, ref);

  epi = linked_ele(ele_info, rt->to_edge);

  if (rt->to_edge->type == 'p' || rt->to_edge->type == 'S') {
    if (((float) ele_info->ele->frag.rb - ele_info->ele->frag.lb)/(epi->ele->frag.rb - epi->ele->frag.lb) < CUTOFF3) {
      if (rt->to_edge->score > ref->to_edge->score) {
        ref->to_edge->type = 'S';
        ref = rt;
      } else {
        rt->to_edge->type = 'S';
      }
    } else {
      if (rt->to_edge->type == 'p') rt->to_edge->type = 'S';
    }
  }

  return ref;
}


void edge_repair(ELE_INFO_t *ele_info) {
    EDGE_TREE_t *hanger;

    ele_read_in(ele_info, 2);
    if (ele_info->ele->edges) {
      hanger = best_link(ele_info, ele_info->ele->edges);
      if (!hanger) {
	fprintf(log_file, "ele %d can not be repaired.\n", ele_info->index);
      } else if (hanger->to_edge->type == 'S') {
	hanger->to_edge->type == 'P';
	ele_write_out(ele_info, 2);
      }
    }
    ele_cleanup(&ele_info->ele);
}


EDGE_TREE_t *best_link(ELE_INFO_t *ele_info, EDGE_TREE_t *rt) {
  EDGE_TREE_t *res=NULL;

  if (rt->to_edge->type == 'p') res = rt;
  else if (rt->to_edge->type == 'P' || rt->to_edge->type == 'S') {
    if (!res || rt->to_edge->score > res->to_edge->score) res = rt;
  }

  if ((!res || res->to_edge->type != 'p') && rt->l) res = best_link(ele_info, rt->l);
  if ((!res || res->to_edge->type != 'p') && rt->r) res = best_link(ele_info, rt->r);

  return res;
}

