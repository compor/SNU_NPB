//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB SP code. This C        //
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

//---------------------------------------------------------------------
// block-diagonal matrix-vector multiplication
//---------------------------------------------------------------------
void pinvr()
{
  int i, j, k;
  double r1, r2, r3, r4, r5, t1, t2;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3;
#endif // USE_CITERATOR

  if (timeron) timer_start(t_pinvr);
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 1, nz2+1, 1, cit_step_add, RND) {
  /*for (k = 1; k <= nz2; k++) {*/
    FOR_START(j, cit2, 1, ny2+1, 1, cit_step_add, RND) {
    /*for (j = 1; j <= ny2; j++) {*/
      FOR_START(i, cit3, 1, nx2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= nx2; i++) {*/
        r1 = rhs[k][j][i][0];
        r2 = rhs[k][j][i][1];
        r3 = rhs[k][j][i][2];
        r4 = rhs[k][j][i][3];
        r5 = rhs[k][j][i][4];

        t1 = bt * r1;
        t2 = 0.5 * ( r4 + r5 );

        rhs[k][j][i][0] =  bt * ( r4 - r5 );
        rhs[k][j][i][1] = -r3;
        rhs[k][j][i][2] =  r2;
        rhs[k][j][i][3] = -t1 + t2;
        rhs[k][j][i][4] =  t1 + t2;
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 1; k <= nz2; k++) {
    for (j = 1; j <= ny2; j++) {
      for (i = 1; i <= nx2; i++) {
        r1 = rhs[k][j][i][0];
        r2 = rhs[k][j][i][1];
        r3 = rhs[k][j][i][2];
        r4 = rhs[k][j][i][3];
        r5 = rhs[k][j][i][4];

        t1 = bt * r1;
        t2 = 0.5 * ( r4 + r5 );

        rhs[k][j][i][0] =  bt * ( r4 - r5 );
        rhs[k][j][i][1] = -r3;
        rhs[k][j][i][2] =  r2;
        rhs[k][j][i][3] = -t1 + t2;
        rhs[k][j][i][4] =  t1 + t2;
      }
    }
  }
#endif // USE_CITERATOR
  if (timeron) timer_stop(t_pinvr);
}

