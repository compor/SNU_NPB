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

#include "header.h"
#include "timers.h"
#include "adt_citerator.h"

#define USE_CITERATOR

//---------------------------------------------------------------------
// addition of update to the vector u
//---------------------------------------------------------------------
void add()
{
  int i, j, k, m;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif

  if (timeron) timer_start(t_add);
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 1, grid_points[2]-2+1, 1, cit_step_add, RND) {
  /*for (k = 1; k <= grid_points[2]-2; k++) {*/
    FOR_START(j, cit2, 1, grid_points[1]-2+1, 1, cit_step_add, RND) {
    /*for (j = 1; j <= grid_points[1]-2; j++) {*/
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        FOR_START(m, cit4, 0, 4+1, 1, cit_step_add, RND) {
        /*for (m = 0; m < 5; m++) {*/
          u[k][j][i][m] = u[k][j][i][m] + rhs[k][j][i][m];
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 1; k <= grid_points[2]-2; k++) {
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 5; m++) {
          u[k][j][i][m] = u[k][j][i][m] + rhs[k][j][i][m];
        }
      }
    }
  }
#endif // USE_CITERATOR
  if (timeron) timer_stop(t_add);
}
