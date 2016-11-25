/* standard stuff */

#include "bolts.h"
#include "seqlist.h"

/* A collection of structures and functions to handle MSPs */

/* BASIC structures */

typedef struct frag {
  char *seq_name;
  int32_t lb, rb;
/*  char lst, rst;    type of ends, can be used to mark physical ends */
} FRAG_t;

typedef struct frag_list {
  FRAG_t *to_frag;
  struct frag_list *next;
} FRAG_DATA_t;

typedef struct image {
  int index;
  FRAG_t frag;
  struct msp *to_msp;
  struct ele_info *ele_info;
} IMAGE_t;

typedef struct image_list {
  IMAGE_t *to_image;
  struct image_list *next;
} IMG_DATA_t;

typedef struct img_tree {
  IMAGE_t *to_image;
  struct img_tree *p, *l, *r;
} IMG_TREE_t;

/*typedef struct img_rec {
  int msp_index;
  struct img_rec *next;
} IMG_REC_t; */

typedef struct msp {
  /*struct msp_list *hanger;*/
  char stat; /* 'p', 's' */
  int32_t score;
  float iden;
  int direction;
  IMAGE_t query, sbjct;
} MSP_t;

typedef struct msp_list {
  MSP_t *to_msp;
  struct msp_list *next;
} MSP_DATA_t;







/* funtions prototypes */

int print_msp(MSP_t *);
int fprint_msp(FILE *, MSP_t *);
int scan_msp(MSP_t *, char *);
int merge_msp(MSP_t *, MSP_t *);
IMAGE_t *partner(IMAGE_t *);
int doub_cov(FRAG_t *, FRAG_t *, float);
int sing_cov(FRAG_t *, FRAG_t *, float);

int print_msp(MSP_t *m) {
  if (m->direction == 1) {
    if (printf("%06d %3.1f %08d %08d %s %08d %08d %s \n", m->score, m->iden, m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name, m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0) {
      return 0;
    } else {
      return 1;
    }
  } else if (m->direction == -1) {
    if (printf("%06d %3.1f %08d %08d %s %08d %08d %s \n", m->score, m->iden, m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name, m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0) {
      return 0;
    } else {
      return 1;
    }
  } else {
    return 1;
  }
}


int fprint_msp(FILE *ofp, MSP_t *m) {
  if (m->direction == 1) {
    if (fprintf(ofp, "%06d %3.1f %08d %08d %s %08d %08d %s \n", m->score, m->iden, m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name, m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0) {
      return 0;
    } else {
      return 1;
    }
  } else if (m->direction == -1) {
    if (fprintf(ofp, "%06d %3.1f %08d %08d %s %08d %08d %s \n", m->score, m->iden, m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name, m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name) > 0) {
      return 0;
    } else {
      return 1;
    }
  } else {
    return 1;
  }
}


int scan_msp(MSP_t *m, char *line) {
    int32_t bd_tmp;
    char qname[NAME_LEN], sname[NAME_LEN];
    int pos;

    if (sscanf(line, "%ld %f %ld %ld %s %ld %ld %s \n", &(m->score), &(m->iden), &(m->query.frag.lb), &(m->query.frag.rb), qname, &(m->sbjct.frag.lb), &(m->sbjct.frag.rb), sname) != 8) {
      return 1;
    } else {
      pos = GetSeqIndex(0, seq_no-1, qname);
      if (pos<0) {return pos;}
      m->query.frag.seq_name = *(seq_names+pos);
      pos = GetSeqIndex(0, seq_no-1, sname);
      if (pos<0) {return pos;}
      m->sbjct.frag.seq_name = *(seq_names+pos);

      m->direction = 1;
      m->stat = 's';
      /*m->score = m->score;*/
      m->query.to_msp = m;
      m->sbjct.to_msp = m;
      /* indices and ele_info pointer are not taken care of here */
      /* m->query.frag.lst = 'n';
      m->query.frag.rst = 'n';
      m->sbjct.frag.lst = 'n';
      m->sbjct.frag.rst = 'n'; */
      if (m->query.frag.lb > m->query.frag.rb) {
	m->direction *= -1;
	bd_tmp = m->query.frag.lb;
	m->query.frag.lb = m->query.frag.rb;
	m->query.frag.rb = bd_tmp;
      }
      if (m->sbjct.frag.lb > m->sbjct.frag.rb) {
	m->direction *= -1;
	bd_tmp = m->sbjct.frag.lb;
	m->sbjct.frag.lb = m->sbjct.frag.rb;
	m->sbjct.frag.rb = bd_tmp;
      }
      return 0;
    }
}

/* archaine */
#if 0
int merge_msp(MSP_t *m1, MSP_t *m2) { /* merge two MSPs into the first one */
  int32_t l1, l2;

  /* no need for changing direction */
  /* iden */
  l1 = m1->sbjct.frag.rb - m1->sbjct.frag.lb;
  l2 = m2->sbjct.frag.rb - m2->sbjct.frag.lb;
  m1->iden = (m1->iden*l1 + m2->iden*l2)/(l1+l2);
  /* query */
  /* links and seq_name do not need to change */
  /* query boundaries */
  m1->query.frag.lb = m2->query.frag.lb < m1->query.frag.lb ? m2->query.frag.lb : m1->query.frag.lb;
  m1->query.frag.rb = m2->query.frag.rb > m1->query.frag.rb ? m2->query.frag.rb : m1->query.frag.rb;
  /* sbjct */
  /* links and seq_name do not need to change */
  /* sbjct boundaries */
  m1->sbjct.frag.lb = m2->sbjct.frag.lb < m1->sbjct.frag.lb ? m2->sbjct.frag.lb : m1->sbjct.frag.lb;
  m1->sbjct.frag.rb = m2->sbjct.frag.rb > m1->sbjct.frag.rb ? m2->sbjct.frag.rb : m1->sbjct.frag.rb;
  /* score */
  if (m1->score == 999999 || m2->score == 999999) {
    m1->score = 999999;
  } else {
    m1->score = (m1->score%800000 + m2->score%800000)*(m1->sbjct.frag.rb - m1->sbjct.frag.lb)/(l1+l2);  /* notice that the boundaries of m1 are now for the merged MSP */
    if (m1->score + 800000 < 999999) {
      m1->score += 800000;
    } else {
      m1->score = 999999;
    }
  }

  return 0;
}
#endif

IMAGE_t *partner(IMAGE_t *i) {
  if (i == &(i->to_msp->query)) {
    return &(i->to_msp->sbjct);
  }
  return &(i->to_msp->query);
}


int sing_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
  int32_t l1, l2, l, lb, rb;
  if (f1->seq_name == f2->seq_name /*!strncmp(f1->seq_name, f2->seq_name, NAME_LEN)*/) {
    l1 = f1->rb - f1->lb;
    l2 = f2->rb - f2->lb;
    lb = f1->lb > f2->lb ? f1->lb : f2->lb;
    rb = f1->rb < f2->rb ? f1->rb : f2->rb;
    l = rb - lb;
    if ( (float) l/l1 >= cutoff || (float) l/l2 >= cutoff) {
      return 1;
    }
  }
  return 0;
}


int doub_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
  int32_t l1, l2, l, lb, rb;
  if (f1->seq_name == f2->seq_name /*!strncmp(f1->seq_name, f2->seq_name, NAME_LEN)*/) {
    l1 = f1->rb - f1->lb;
    l2 = f2->rb - f2->lb;
    lb = f1->lb > f2->lb ? f1->lb : f2->lb;
    rb = f1->rb < f2->rb ? f1->rb : f2->rb;
    l = rb - lb;
    if ( (float) l/l1 >= cutoff && (float) l/l2 >= cutoff) {
      return 1;
    }
  }
  return 0;
}
