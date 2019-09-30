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

/****/

#define FOR_START(IT, CIT, LB, UB, STEP, ORDER) \
    cit_create(&(CIT), (LB), (UB), (STEP), (ORDER)); \
    for (IT = cit_begin((CIT)); cit_is_valid((CIT)); IT = cit_next((CIT)))

#define FOR_END(CIT) \
    do { cit_destroy((CIT)); } while(0)

#define FOR_RND_START(IT, CIT, LB, UB, STEP) \
    FOR_START((IT), (CIT), (LB), (UB), (STEP), (RND)) \

#define FOR_RND_END(CIT) \
    FOR_END((CIT))

#endif /* ADT_CITERATOR_H */
