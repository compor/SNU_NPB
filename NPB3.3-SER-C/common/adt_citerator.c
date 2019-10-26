
#include "adt_citerator.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

struct cit_data {
  cit_int_t start;
  cit_int_t end;
  cit_int_t current_idx;
  cit_int_t remaining_indices;

  short int valid;

  cit_int_t *indices;
  cit_int_t nindices;

  int (*stepfunc)(const struct cit_data *data);
  int step;

  enum cit_order order;
};

/****/

void cit_create(struct cit_data **data, cit_int_t start, cit_int_t end, int step,
  int (*stepfunc)(const struct cit_data *data), enum cit_order order) {
  struct timeval time;

  *data = (struct cit_data *) malloc(1 * sizeof(struct cit_data));
  if(!*data) {
    /* TODO add message */
    return;
  }

  (*data)->stepfunc = stepfunc;
  (*data)->step= step;

  if(start < end) {
    (*data)->nindices = (end - start) / abs((*data)->stepfunc(*data));

    if((end - start) % abs((*data)->stepfunc(*data))) {
      (*data)->nindices++;
    }
  }
  else if (start > end) {
    (*data)->nindices = (start - end) / abs((*data)->stepfunc(*data));

    if((start - end) % abs((*data)->stepfunc(*data))) {
      (*data)->nindices++;
    }
  }
  else {
    (*data)->nindices = 0;
  }

  (*data)->remaining_indices = (*data)->nindices;

  if(RND == order && (*data)->nindices) {
    (*data)->indices = malloc((*data)->nindices * sizeof(cit_int_t));
    if(!(*data)->indices) {
      /* TODO add message */
      return;
    }

    (*data)->indices[0] = start;
    cit_int_t i = 0;
    for(i = 1; i < (*data)->nindices; i++) {
      (*data)->indices[i] = (*data)->indices[i-1] + (*data)->stepfunc(*data);
    }

    gettimeofday(&time,NULL);
    srand((time.tv_sec * 1000) + (time.tv_usec / 1000));
    /*srand(time(NULL));*/

    /* randomly permute indices */
    cit_int_t j = 0;
    cit_int_t temp = 0;
    for (i = (*data)->nindices-1; i > 0; --i){
      //generate a random number [0, n-1]
      j = (cit_int_t) rand() % (i+1);

      //swap the last element with element at random index
      temp = (*data)->indices[i];
      (*data)->indices[i] = (*data)->indices[j];
      (*data)->indices[j] = temp;
    }

    (*data)->start = end;
    (*data)->end = start;
    // repurpose
    (*data)->current_idx = 0;
  } else if(BWD == order) {
    (*data)->start = end;
    (*data)->end = start;
    (*data)->current_idx = (*data)->start;
  }
  else {
    (*data)->start = start;
    (*data)->end = end;
    (*data)->current_idx = (*data)->start;
    order = FWD;
  }

  (*data)->valid = ((*data)->remaining_indices) ? 1 : 0;
  (*data)->order = order;
}

void cit_destroy(struct cit_data *data) {
  if(data) {
    if(RND == data->order && data->indices) {
      free(data->indices);
    }
    free(data);
  }
}

cit_int_t cit_begin(struct cit_data *data) {
  if(RND == data->order && data->nindices) {
    return data->indices[data->current_idx];
  }
  else {
   return data->start;
  }
}

cit_int_t cit_next(struct cit_data *data) {
  int stepfunc = data->stepfunc(data);

  if(data->remaining_indices <= 1) {
    data->remaining_indices = 0;
    data->valid = 0;

    if(RND == data->order && data->nindices) {
      return data->indices[data->nindices - 1];
    }
    else {
      return data->end;
    }
  }

  --data->remaining_indices;

  if(RND == data->order) {
    return data->indices[++data->current_idx];
  }
  else if(BWD == data->order) {
    data->current_idx -= stepfunc;
  }
  else {
    data->current_idx += stepfunc;
  }

  return data->current_idx;
}

short cit_is_valid(const struct cit_data *data) {
  return data->valid;
}

/****/

int cit_step_add(const struct cit_data *data) {
  return data->step;
}

