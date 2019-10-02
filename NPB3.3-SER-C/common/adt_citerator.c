
#include "adt_citerator.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

struct cit_data {
  unsigned start;
  unsigned end;
  unsigned current_idx;
  unsigned remaining_indices;

  short int valid;

  unsigned *indices;
  unsigned nindices;

  int (*step)(const struct cit_data *data);

  enum cit_order order;
};

/****/

void cit_create(struct cit_data **data, unsigned start, unsigned end,
  int (*step)(const struct cit_data *data), enum cit_order order) {

  *data = (struct cit_data *) malloc(1 * sizeof(struct cit_data));
  if(!*data) {
    /* TODO add message */
    return;
  }

  (*data)->step = step;

  if(start <= end) {
    (*data)->nindices = (end - start) / abs((*data)->step(*data));
    if((end - start) / abs((*data)->step(*data))) {
      (*data)->nindices++;
    }
  }
  else {
    (*data)->nindices = (start - end) / abs((*data)->step(*data));
    if((start - end) / abs((*data)->step(*data))) {
      (*data)->nindices++;
    }
  }

  (*data)->remaining_indices = (*data)->nindices;

  if(RND == order) {
    (*data)->indices = malloc((*data)->nindices * sizeof(unsigned));
    if(!(*data)->indices) {
      /* TODO add message */
      return;
    }

    (*data)->indices[0] = start;
    unsigned i = 0;
    for(i = 1; i < (*data)->nindices; i++) {
      (*data)->indices[i] = (*data)->indices[i-1] + (*data)->step(*data);
    }

    srand(time(NULL));

    /* randomly permute indices */
    unsigned j = 0;
    unsigned temp = 0;
    for (i = (*data)->nindices-1; i > 0; --i){
      //generate a random number [0, n-1]
      j = (unsigned) rand() % (i+1);

      //swap the last element with element at random index
      temp = (*data)->indices[i];
      (*data)->indices[i] = (*data)->indices[j];
      (*data)->indices[j] = temp;
    }

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

unsigned cit_begin(struct cit_data *data) {
  if(RND == data->order) {
    return data->indices[data->current_idx];
  }
  else {
   return data->start;
  }
}

unsigned cit_next(struct cit_data *data) {
  int step = data->step(data);

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
    data->current_idx -= step;
  }
  else {
    data->current_idx += step;
  }

  return data->current_idx;
}

short cit_is_valid(const struct cit_data *data) {
  return data->valid;
}

/****/

int cit_inc1(const struct cit_data *data) {
  return 1;
}

int cit_dec1(const struct cit_data *data) {
  return -1;
}

