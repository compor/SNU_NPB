//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB UA code. This C        //
//  version is developed by the Center for Manycore Programming at Seoul   //
//  National University and derived from the serial Fortran versions in    //
//  "NPB3.3-SER" developed by NAS.                                         //
//                                                                         //
//  Permission to use, copy, distribute and modify this software for any   //
//  purpose with or without fee is hereby granted. This software is        //
//  provided "as is" without express or implied warranty.                  //
//                                                                         //
//  Information on NPB 3.3, including the technical report, the original   //
//  specifications, source code, results and information on how to submit  //
//  new results, is available at:                                          //
//                                                                         //
//           http://www.nas.nasa.gov/Software/NPB/                         //
//                                                                         //
//  Send comments or suggestions for this C version to cmp@aces.snu.ac.kr  //
//                                                                         //
//          Center for Manycore Programming                                //
//          School of Computer Science and Engineering                     //
//          Seoul National University                                      //
//          Seoul 151-744, Korea                                           //
//                                                                         //
//          E-mail:  cmp@aces.snu.ac.kr                                    //
//                                                                         //
//-------------------------------------------------------------------------//

//-------------------------------------------------------------------------//
// Authors: Sangmin Seo, Jungwon Kim, Jun Lee, Jeongho Nah, Gangwon Jo,    //
//          and Jaejin Lee                                                 //
//-------------------------------------------------------------------------//

#include "header.h"
#include "adt_citerator.h"

#define USE_CITERATOR

//------------------------------------------------------------------
// initialize double precision array a with length of n
//------------------------------------------------------------------
void reciprocal(double a[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = 1.0/a[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = 1.0/a[i];
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// initialize double precision array a with length of n
//------------------------------------------------------------------
void r_init(double a[], int n, double _const)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = _const;
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = _const;
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// initialize integer array a with length of n
//------------------------------------------------------------------
void nr_init(int a[], int n, int _const)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = _const;
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = _const;
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// initialize logical array a with length of n
//------------------------------------------------------------------
void l_init(logical a[], int n, logical _const)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = _const;
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = _const;
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// copy array of integers b to a, the length of array is n
//------------------------------------------------------------------
void ncopy(int a[], int b[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = b[i];
  }
#endif // USE_CITERATOR

}


//------------------------------------------------------------------
// copy double precision array b to a, the length of array is n
//------------------------------------------------------------------
void copy(double a[], double b[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = b[i];
  }
#endif // USE_CITERATOR
}


//-----------------------------------------------------------------
// a=b*c1
//-----------------------------------------------------------------
void adds2m1(double a[], double b[], double c1, int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = a[i]+c1*b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = a[i]+c1*b[i];
  }
#endif // USE_CITERATOR
}


//-----------------------------------------------------------------
// a=c1*a+b
//-----------------------------------------------------------------
void adds1m1(double a[], double b[], double c1, int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = c1*a[i]+b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = c1*a[i]+b[i];
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// a=a*b
//------------------------------------------------------------------
void col2(double a[], double b[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = a[i]*b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = a[i]*b[i];
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// zero out array of integers
//------------------------------------------------------------------
void nrzero(int na[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    na[i] = 0;
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    na[i] = 0;
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// a=a+b
//------------------------------------------------------------------
void add2(double a[], double b[], int n)
{
  int i;
#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(i, cit1, 0, n, 1, cit_step_add, RND) {
  /*for (i = 0; i < n; i++) {*/
    a[i] = a[i]+b[i];
  }
  FOR_END(cit1);
#else
  for (i = 0; i < n; i++) {
    a[i] = a[i]+b[i];
  }
#endif // USE_CITERATOR
}


//------------------------------------------------------------------
// calculate the integral of ta1 over the whole domain
//------------------------------------------------------------------
double calc_norm()
{
  double total, ieltotal;
  int iel, k, j, i, isize;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif // USE_CITERATOR

  total = 0.0;

#ifdef USE_CITERATOR
  FOR_START(iel, cit1, 0, nelt, 1, cit_step_add, RND) {
  /*for (iel = 0; iel < nelt; iel++) {*/
    ieltotal = 0.0;
    isize = size_e[iel];
    FOR_START(k, cit2, 0, LX1, 1, cit_step_add, RND) {
    /*for (k = 0; k < LX1; k++) {*/
      FOR_START(j, cit3, 0, LX1, 1, cit_step_add, RND) {
      /*for (j = 0; j < LX1; j++) {*/
        FOR_START(i, cit4, 0, LX1, 1, cit_step_add, RND) {
        /*for (i = 0; i < LX1; i++) {*/
          ieltotal = ieltotal+ta1[iel][k][j][i]*w3m1[k][j][i]
                    *jacm1_s[isize][k][j][i];
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
    total = total+ieltotal;
  }
  FOR_END(cit1);
#else
  for (iel = 0; iel < nelt; iel++) {
    ieltotal = 0.0;
    isize = size_e[iel];
    for (k = 0; k < LX1; k++) {
      for (j = 0; j < LX1; j++) {
        for (i = 0; i < LX1; i++) {
          ieltotal = ieltotal+ta1[iel][k][j][i]*w3m1[k][j][i]
                    *jacm1_s[isize][k][j][i];
        }
      }
    }
    total = total+ieltotal;
  }
#endif // USE_CITERATOR

  return total;
}


//-----------------------------------------------------------------
// input array frontier, perform (potentially) parallel add so that
// the output frontier[i] has sum of frontier[1]+frontier[2]+...+frontier[i]
//-----------------------------------------------------------------
void parallel_add(int frontier[])
{
  int nellog, i, ahead, ii, ntemp, n1, ntemp1, n2, iel;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3;
#endif // USE_CITERATOR

  if (nelt <= 1) return;

  nellog = 0;
  iel = 1;
  do {
    iel = iel*2;
    nellog = nellog+1;
  } while (iel < nelt);

  ntemp = 1;
#ifdef _USE_CITERATOR
  FOR_START(i, cit1, 0, nellog, 1, cit_step_add, RND) {
  /*for (i = 0; i < nellog; i++) {*/
    /*fprintf(stderr, "lol\n");*/
    n1 = ntemp*2;
    n2 = n1;
    FOR_START(iel, cit2, n1, nelt+1, n1, cit_step_add, RND) {
    /*for (iel = n1; iel <= nelt; iel += n1) {*/
      /*fprintf(stderr, "lal\n");*/
      ahead = frontier[iel-ntemp-1];
      FOR_START(ii, cit3, ntemp-1, -1, -1, cit_step_add, RND) {
      /*for (ii = ntemp-1; ii >= 0; ii--) {*/
        frontier[iel-ii-1] = frontier[iel-ii-1]+ahead;
      }
      FOR_END(cit3);
      n2 = iel;
    }
    FOR_END(cit2);
    if (n2 <= nelt) n2 = n2+n1;

    ntemp1 = n2-nelt;
    if (ntemp1 < ntemp) {
      ahead = frontier[n2-ntemp-1];
      FOR_START(ii, cit3, ntemp-1, ntemp1-1, -1, cit_step_add, RND) {
      /*for (ii = ntemp-1; ii >= ntemp1; ii--) {*/
        frontier[n2-ii-1] = frontier[n2-ii-1]+ahead;
      }
      FOR_END(cit3);
    }

    ntemp = n1;
  }
  FOR_END(cit1);
#else
  for (i = 0; i < nellog; i++) {
    n1 = ntemp*2;
    n2 = n1;
    for (iel = n1; iel <= nelt; iel += n1) {
      ahead = frontier[iel-ntemp-1];
      for (ii = ntemp-1; ii >= 0; ii--) {
        frontier[iel-ii-1] = frontier[iel-ii-1]+ahead;
      }
      n2 = iel;
    }
    if (n2 <= nelt) n2 = n2+n1;

    ntemp1 = n2-nelt;
    if (ntemp1 < ntemp) {
      ahead = frontier[n2-ntemp-1];
      for (ii = ntemp-1; ii >= ntemp1; ii--) {
        frontier[n2-ii-1] = frontier[n2-ii-1]+ahead;
      }
    }

    ntemp = n1;
  }
#endif // USE_CITERATOR
}



//------------------------------------------------------------------
// Perform stiffness summation: element-mortar-element mapping
//------------------------------------------------------------------
void dssum()
{
  transfb(dpcmor, (double *)dpcelm);
  transf (dpcmor, (double *)dpcelm);
}


//------------------------------------------------------------------
// assign the value val to face(iface,iel) of array a.
//------------------------------------------------------------------
void facev(double a[LX1][LX1][LX1], int iface, double val)
{
  int kx1, kx2, ky1, ky2, kz1, kz2, ix, iy, iz;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3;
#endif // USE_CITERATOR

  kx1 = 1;
  ky1 = 1;
  kz1 = 1;
  kx2 = LX1;
  ky2 = LX1;
  kz2 = LX1;
  if (iface == 0) kx1 = LX1;
  if (iface == 1) kx2 = 1;
  if (iface == 2) ky1 = LX1;
  if (iface == 3) ky2 = 1;
  if (iface == 4) kz1 = LX1;
  if (iface == 5) kz2 = 1;

#ifdef USE_CITERATOR
  FOR_START(ix, cit1, kx1-1, kx2, 1, cit_step_add, RND) {
  /*for (ix = kx1-1; ix < kx2; ix++) {*/
    FOR_START(iy, cit2, ky1-1, ky2, 1, cit_step_add, RND) {
    /*for (iy = ky1-1; iy < ky2; iy++) {*/
      FOR_START(iz, cit3, kz1-1, kz2, 1, cit_step_add, RND) {
      /*for (iz = kz1-1; iz < kz2; iz++) {*/
        a[iz][iy][ix] = val;
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (ix = kx1-1; ix < kx2; ix++) {
    for (iy = ky1-1; iy < ky2; iy++) {
      for (iz = kz1-1; iz < kz2; iz++) {
        a[iz][iy][ix] = val;
      }
    }
  }
#endif // USE_CITERATOR
}
