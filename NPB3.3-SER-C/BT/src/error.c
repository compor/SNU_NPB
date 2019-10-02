//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB BT code. This C        //
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

#include <math.h>
#include "header.h"
#include "adt_citerator.h"

#define USE_CITERATOR

//---------------------------------------------------------------------
// this function computes the norm of the difference between the
// computed solution and the exact solution
//---------------------------------------------------------------------
void error_norm(double rms[5])
{
  int i, j, k, m, d;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif // USE_CITERATOR
  double xi, eta, zeta, u_exact[5], add;

#ifdef USE_CITERATOR
  FOR_RND_START(m, cit1, 0, 4, CIT_STEP1) {
  /*for (m = 0; m < 5; m++) {*/
    rms[m] = 0.0;
  }
  FOR_RND_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(k, cit1, 0, grid_points[2]-1, CIT_STEP1) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_RND_START(j, cit2, 0, grid_points[1]-1, CIT_STEP1) {
    /*for (j = 0; j <= grid_points[1]-1; j++) {*/
      eta = (double)(j) * dnym1;
      FOR_RND_START(i, cit3, 0, grid_points[0]-1, CIT_STEP1) {
      /*for (i = 0; i <= grid_points[0]-1; i++) {*/
        xi = (double)(i) * dnxm1;
        exact_solution(xi, eta, zeta, u_exact);

        FOR_RND_START(m, cit4, 0, 4, CIT_STEP1) {
        /*for (m = 0; m < 5; m++) {*/
          add = u[k][j][i][m]-u_exact[m];
          rms[m] = rms[m] + add*add;
        }
        FOR_RND_END(cit4);
      }
      FOR_RND_END(cit3);
    }
    FOR_RND_END(cit2);
  }
  FOR_RND_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      for (i = 0; i <= grid_points[0]-1; i++) {
        xi = (double)(i) * dnxm1;
        exact_solution(xi, eta, zeta, u_exact);

        for (m = 0; m < 5; m++) {
          add = u[k][j][i][m]-u_exact[m];
          rms[m] = rms[m] + add*add;
        }
      }
    }
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(m, cit1, 0, 4, CIT_STEP1) {
  /*for (m = 0; m < 5; m++) {*/
    FOR_RND_START(d, cit2, 0, 2, CIT_STEP1) {
    /*for (d = 0; d < 3; d++) {*/
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    FOR_RND_END(cit2);
    rms[m] = sqrt(rms[m]);
  }
  FOR_RND_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
#endif // USE_CITERATOR
}


void rhs_norm(double rms[5])
{
  int i, j, k, d, m;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif // USE_CITERATOR
  double add;

#ifdef USE_CITERATOR
  FOR_RND_START(m, cit1, 0, 4, CIT_STEP1) {
  /*for (m = 0; m < 5; m++) {*/
    rms[m] = 0.0;
  }
  FOR_RND_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(k, cit1, 1, grid_points[2]-2, CIT_STEP1) {
  /*for (k = 1; k <= grid_points[2]-2; k++) {*/
    FOR_RND_START(j, cit2, 1, grid_points[1]-2, CIT_STEP1) {
    /*for (j = 1; j <= grid_points[1]-2; j++) {*/
      FOR_RND_START(i, cit3, 1, grid_points[0]-2, CIT_STEP1) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        FOR_RND_START(m, cit4, 0, 4, CIT_STEP1) {
        /*for (m = 0; m < 5; m++) {*/
          add = rhs[k][j][i][m];
          rms[m] = rms[m] + add*add;
        }
        FOR_RND_END(cit4);
      }
      FOR_RND_END(cit3);
    }
    FOR_RND_END(cit2);
  }
  FOR_RND_END(cit1);
#else
  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 5; m++) {
          add = rhs[k][j][i][m];
          rms[m] = rms[m] + add*add;
        }
      }
    }
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_RND_START(m, cit1, 0, 4, CIT_STEP1) {
  /*for (m = 0; m < 5; m++) {*/
    FOR_RND_START(d, cit2, 0, 2, CIT_STEP1) {
    /*for (d = 0; d < 3; d++) {*/
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    FOR_RND_END(cit2);
    rms[m] = sqrt(rms[m]);
  }
  FOR_RND_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
#endif // USE_CITERATOR
}
