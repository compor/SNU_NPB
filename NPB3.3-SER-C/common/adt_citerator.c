
#include "adt_citerator.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

struct cit_data {
  unsigned int start;
  unsigned int end;
  unsigned int current;

  short int valid;

  unsigned int *indices;
  unsigned int nindices;

  unsigned int (*step)(const struct cit_data *data);

  enum cit_order order;
};

/****/

void cit_create(struct cit_data **data, unsigned int start, unsigned int end,
  unsigned (*step)(const struct cit_data *data), enum cit_order order) {

  *data = (struct cit_data *) malloc(1 * sizeof(struct cit_data));
  if(!*data) {
    /* TODO add message */
    return;
  }

  if(RND == order) {
    (*data)->start = start;
    (*data)->end = end;
    (*data)->step = step;

    (*data)->nindices = (end - start) / (*data)->step(*data);
    if((end - start) / (*data)->step(*data)) {
      (*data)->nindices++;
    }

    (*data)->indices = malloc((*data)->nindices * sizeof(unsigned));
    /* TODO missing check */

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
    (*data)->current = 0;
  } else if(BWD == order) {
    (*data)->start = end;
    (*data)->end = start;
    (*data)->step = step;
    (*data)->current = (*data)->start;
  }
  else {
    (*data)->start = start;
    (*data)->end = end;
    (*data)->step = step;
    (*data)->current = (*data)->start;
    order = FWD;
  }

  (*data)->valid = 1;
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

unsigned int cit_begin(struct cit_data *data) {
  if(RND == data->order) {
    return data->indices[data->current];
  }
  else {
   return data->start;
  }
}

unsigned int cit_next(struct cit_data *data) {
  unsigned int step = data->step(data);

  if(RND == data->order) {
    return data->indices[++data->current];
  } else if(BWD == data->order) {
    if(data->current == data->end) {
      data->valid = 0;
    } else {
      data->current -= step;
    }
  }
  else {
    data->current += step;
  }

  return data->current;
}

short cit_is_valid(const struct cit_data *data) {
  if(RND == data->order) {
    return data->current < data->nindices;
  } else if(BWD == data->order) {
    return data->valid && data->current >= data->end;
  }
  else {
    return data->current <= data->end;
  }
}

/****/

unsigned int cit_step1(const struct cit_data *data) {
  return 1;
}

