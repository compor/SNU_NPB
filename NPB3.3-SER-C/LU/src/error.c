//-------------------------------------------------------------------------//
//                                                                         //
//  This benchmark is a serial C version of the NPB LU code. This C        //
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
#include "applu.incl"
#include "adt_citerator.h"

#define USE_CITERATOR

//---------------------------------------------------------------------
//
// compute the solution error
//
//---------------------------------------------------------------------
void error()
{
  //---------------------------------------------------------------------
  // local variables
  //---------------------------------------------------------------------
  int i, j, k, m;
  double tmp;
  double u000ijk[5];
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(m, cit1, 0, 5, 1, cit_step_add, RND) {
  /*for (m = 0; m < 5; m++) {*/
    errnm[m] = 0.0;
  }
  FOR_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    errnm[m] = 0.0;
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(k, cit1, 1, nz-1, 1, cit_step_add, RND) {
  /*for (k = 1; k < nz-1; k++) {*/
    FOR_START(j, cit2, jst, jend, 1, cit_step_add, RND) {
    /*for (j = jst; j < jend; j++) {*/
      FOR_START(i, cit3, ist, iend, 1, cit_step_add, RND) {
      /*for (i = ist; i < iend; i++) {*/
        exact( i, j, k, u000ijk );
        FOR_START(m, cit4, 0, 5, 1, cit_step_add, RND) {
        /*for (m = 0; m < 5; m++) {*/
          tmp = ( u000ijk[m] - u[k][j][i][m] );
          errnm[m] = errnm[m] + tmp * tmp;
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 1; k < nz-1; k++) {
    for (j = jst; j < jend; j++) {
      for (i = ist; i < iend; i++) {
        exact( i, j, k, u000ijk );
        for (m = 0; m < 5; m++) {
          tmp = ( u000ijk[m] - u[k][j][i][m] );
          errnm[m] = errnm[m] + tmp * tmp;
        }
      }
    }
  }
#endif // USE_CITERATOR

#ifdef USE_CITERATOR
  FOR_START(m, cit1, 0, 5, 1, cit_step_add, RND) {
  /*for (m = 0; m < 5; m++) {*/
    errnm[m] = sqrt ( errnm[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
  FOR_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    errnm[m] = sqrt ( errnm[m] / ( (nx0-2)*(ny0-2)*(nz0-2) ) );
  }
#endif // USE_CITERATOR

  /*
  printf(" \n RMS-norm of error in soln. to first pde  = %12.5E\n"
         " RMS-norm of error in soln. to second pde = %12.5E\n"
         " RMS-norm of error in soln. to third pde  = %12.5E\n"
         " RMS-norm of error in soln. to fourth pde = %12.5E\n"
         " RMS-norm of error in soln. to fifth pde  = %12.5E\n",
         errnm[0], errnm[1], errnm[2], errnm[3], errnm[4]);
  */
}

