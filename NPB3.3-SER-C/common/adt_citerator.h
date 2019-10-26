#ifndef ADT_CITERATOR_H
#define ADT_CITERATOR_H

struct cit_data;

enum cit_order {
  FWD,
  BWD,
  RND, /* NOT IMPLEMENTED */
};

typedef int cit_int_t;

/****/

void cit_create(struct cit_data **data, cit_int_t start, cit_int_t end,
                int step, int (*stepfunc)(const struct cit_data *data),
                enum cit_order order);
void cit_destroy(struct cit_data *data);

cit_int_t cit_begin(struct cit_data *data);
cit_int_t cit_next(struct cit_data *data);
short cit_is_valid(const struct cit_data *data);

/****/

int cit_step_add(const struct cit_data *data);

/****/

#define FOR_START(IT, CIT, LB, UB, STEP, STEPVAL, ORDER) \
    cit_create(&(CIT), (LB), (UB), (STEP), (STEPVAL), (ORDER)); \
    for (IT = cit_begin((CIT)); cit_is_valid((CIT)); IT = cit_next((CIT)))

#define FOR_END(CIT) \
    do { cit_destroy((CIT)); } while(0)

#define FOR_RND_START(IT, CIT, LB, UB, STEP, STEPVAL) \
    FOR_START((IT), (CIT), (LB), (UB), (STEP), (STEPVAL), (RND)) \

#define FOR_RND_END(CIT) \
    FOR_END((CIT))

#endif /* ADT_CITERATOR_H */
