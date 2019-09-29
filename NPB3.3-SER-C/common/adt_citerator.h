#ifndef ADT_CITERATOR_H
#define ADT_CITERATOR_H

struct cit_data;

enum cit_order {
  FWD,
  BWD,
  RND, /* NOT IMPLEMENTED */
};

/****/

void cit_create(struct cit_data **data, unsigned int start, unsigned int end,
                unsigned (*step)(const struct cit_data *data),
                enum cit_order order);
void cit_destroy(struct cit_data *data);

unsigned int cit_begin(struct cit_data *data);
unsigned int cit_next(struct cit_data *data);
short cit_is_valid(const struct cit_data *data);

/****/

unsigned int cit_step1(const struct cit_data *data);

#define CIT_STEP1 cit_step1

#endif /* ADT_CITERATOR_H */
