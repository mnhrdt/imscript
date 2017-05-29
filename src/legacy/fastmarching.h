
#ifndef _FASTMARCHING_H
#define _FASTMARCHING_H

//#include "qnm.h"
#include "imatge3d.h"

void eikonal_slow(imatge *out, imatge *in);
void fill_distance_fast(imatge *out, float (*d)[3], int n);
void fill_distance_extemely_slow(imatge *out, float (*d)[3], int n);
void build_signed_distance(imatge *, imatge *);
void fill_hinted_signed_distance(imatge *, float (*)[3], int, imatge *);

#endif /* _FASTMARCHING_H */
