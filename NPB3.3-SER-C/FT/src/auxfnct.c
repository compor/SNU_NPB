//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB FT code. This C        //
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

#include <stdio.h>
#include <math.h>

#include "global.h"
#include "randdp.h"

#include "adt_citerator.h"

#define USE_CITERATOR

//---------------------------------------------------------------------
// compute the roots-of-unity array that will be used for subsequent FFTs.
//---------------------------------------------------------------------
void CompExp(int n, dcomplex exponent[n])
{
  int m, nu, ku, i, j, ln;
  double t, ti;
  const double pi = 3.141592653589793238;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2;
#endif // USE_CITERATOR

  nu = n;
  m = ilog2(n);
  exponent[0] = dcmplx(m, 0.0);
  ku = 2;
  ln = 1;
#ifdef USE_CITERATOR
  FOR_RND_START(j, cit1, 1, m, 1, cit_step_add) {
  /*for (j = 1; j <= m; j++) {*/
    t = pi / ln;
    FOR_RND_START(i, cit2, 0, ln - 1, 1, cit_step_add) {
    /*for (i = 0; i <= ln - 1; i++) {*/
      ti = i * t;
      exponent[i+ku-1] = dcmplx(cos(ti), sin(ti));
    }
    FOR_RND_END(cit2);
    ku = ku + ln;
    ln = 2 * ln;
  }
  FOR_RND_END(cit1);
#else
  for (j = 1; j <= m; j++) {
    t = pi / ln;
    for (i = 0; i <= ln - 1; i++) {
      ti = i * t;
      exponent[i+ku-1] = dcmplx(cos(ti), sin(ti));
    }
    ku = ku + ln;
    ln = 2 * ln;
  }
#endif // USE_CITERATOR
}


int ilog2(int n)
{
  int nn, lg;
  if (n == 1) return 0;

  lg = 1;
  nn = 2;
  while (nn < n) {
    nn = nn * 2;
    lg = lg + 1;
  }
  return lg;
}


//---------------------------------------------------------------------
// compute a^exponent mod 2^46
//---------------------------------------------------------------------
static double ipow46(double a, int exponent)
{
  double result, dummy, q, r;
  int n, n2;

  //---------------------------------------------------------------------
  // Use
  //   a^n = a^(n/2)*a^(n/2) if n even else
  //   a^n = a*a^(n-1)       if n odd
  //---------------------------------------------------------------------
  result = 1;
  if (exponent == 0) return result;
  q = a;
  r = 1;
  n = exponent;

  while (n > 1) {
    n2 = n / 2;
    if (n2 * 2 == n) {
      dummy = randlc(&q, q);
      n = n2;
    } else {
      dummy = randlc(&r, q);
      n = n-1;
    }
  }
  dummy = randlc(&r, q);
  result = r;
  return result;
}


void CalculateChecksum(dcomplex *csum, int iterN, int d1, int d2, int d3,
                       dcomplex u[d3][d2][d1+1])
{
  int i, i1, ii, ji, ki;

#ifdef USE_CITERATOR
  struct cit_data *cit1;
#endif // USE_CITERATOR

  dcomplex csum_temp = dcmplx(0.0, 0.0);
#ifdef USE_CITERATOR
  FOR_RND_START(i, cit1, 1, 1024, 1, cit_step_add) {
  /*for (i = 1; i <= 1024; i++) {*/
    i1 = i;
    ii = i1 % d1;
    ji = 3 * i1 % d2;
    ki = 5 * i1 % d3;
    csum_temp = dcmplx_add(csum_temp, u[ki][ji][ii]);
  }
  FOR_RND_END(cit1);
#else
  for (i = 1; i <= 1024; i++) {
    i1 = i;
    ii = i1 % d1;
    ji = 3 * i1 % d2;
    ki = 5 * i1 % d3;
    csum_temp = dcmplx_add(csum_temp, u[ki][ji][ii]);
  }
#endif // USE_CITERATOR
  csum_temp = dcmplx_div2(csum_temp, (double)(d1*d2*d3));
  printf(" T =%5d     Checksum =%22.12E%22.12E\n",
      iterN, csum_temp.real, csum_temp.imag);
  *csum = csum_temp;
}


void compute_initial_conditions(int d1, int d2, int d3,
                                dcomplex u0[d3][d2][d1+1])
{
  dcomplex tmp[MAXDIM];
  double x0, start, an, dummy;
  double RanStarts[MAXDIM];

  int i, j, k;

#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3;
#endif // USE_CITERATOR

  const double seed = 314159265.0;
  const double a = 1220703125.0;

  start = seed;
  //---------------------------------------------------------------------
  // Jump to the starting element for our first plane.
  //---------------------------------------------------------------------
  an = ipow46(a, 0);
  dummy = randlc(&start, an);
  an = ipow46(a, 2*d1*d2);
  //---------------------------------------------------------------------
  // Go through by z planes filling in one square at a time.
  //---------------------------------------------------------------------
  RanStarts[0] = start;
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 1, d3-1, 1, cit_step_add, FWD) {
  /*for (k = 1; k < d3; k++) {*/
    dummy = randlc(&start, an);
    RanStarts[k] = start;
  }
  FOR_END(cit1);
#else
  for (k = 1; k < d3; k++) {
    dummy = randlc(&start, an);
    RanStarts[k] = start;
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(k, cit1, 0, d3-1, 1, cit_step_add) {
  /*for (k = 0; k < d3; k++) {*/
    x0 = RanStarts[k];
    FOR_START(j, cit2, 0, d2-1, 1, cit_step_add, FWD) {
    /*for (j = 0; j < d2; j++) {*/
      vranlc(2*d1, &x0, a, (double *)tmp);
      FOR_RND_START(i, cit3, 0, d1-1, 1, cit_step_add) {
      /*for (i = 0; i < d1; i++) {*/
        u0[k][j][i] = tmp[i];
      }
      FOR_RND_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_RND_END(cit1);
#else
  for (k = 0; k < d3; k++) {
    x0 = RanStarts[k];
    for (j = 0; j < d2; j++) {
      vranlc(2*d1, &x0, a, (double *)tmp);
      for (i = 0; i < d1; i++) {
        u0[k][j][i] = tmp[i];
      }
    }
  }
#endif // USE_CITERATOR
}


void evolve(int nx, int ny, int nz,
            dcomplex x[nz][ny][nx+1], dcomplex y[nz][ny][nx+1],
            double twiddle[nz][ny][nx+1])
{
  int i, j, k;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(i, cit1, 0, nz-1, 1, cit_step_add) {
  /*for (i = 0; i < nz; i++) {*/
    FOR_RND_START(k, cit2, 0, ny-1, 1, cit_step_add) {
    /*for (k = 0; k < ny; k++) {*/
      FOR_RND_START(j, cit3, 0, nx-1, 1, cit_step_add) {
      /*for (j = 0; j < nx; j++) {*/
        y[i][k][j] = dcmplx_mul2(y[i][k][j], twiddle[i][k][j]);
        x[i][k][j] = y[i][k][j];
      }
      FOR_RND_END(cit3);
    }
    FOR_RND_END(cit2);
  }
  FOR_RND_END(cit1);
#else
  for (i = 0; i < nz; i++) {
    for (k = 0; k < ny; k++) {
      for (j = 0; j < nx; j++) {
        y[i][k][j] = dcmplx_mul2(y[i][k][j], twiddle[i][k][j]);
        x[i][k][j] = y[i][k][j];
      }
    }
  }
#endif // USE_CITERATOR
}

