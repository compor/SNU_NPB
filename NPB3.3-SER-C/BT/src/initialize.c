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
#include "adt_citerator.h"

#define USE_CITERATOR

//---------------------------------------------------------------------
// This subroutine initializes the field variable u using
// tri-linear transfinite interpolation of the boundary values
//---------------------------------------------------------------------
void initialize()
{
  int i, j, k, m, ix, iy, iz;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

  //---------------------------------------------------------------------
  // Later (in compute_rhs) we compute 1/u for every element. A few of
  // the corner elements are not used, but it convenient (and faster)
  // to compute the whole thing with a simple loop. Make sure those
  // values are nonzero by initializing the whole thing here.
  //---------------------------------------------------------------------
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    FOR_START(j, cit2, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
    /*for (j = 0; j <= grid_points[1]-1; j++) {*/
      FOR_START(i, cit3, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
      /*for (i = 0; i <= grid_points[0]-1; i++) {*/
        FOR_START(m, cit4, 0, 4+1, 1, cit_step_add, RND) {
        /*for (m = 0; m < 5; m++) {*/
          u[k][j][i][m] = 1.0;
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    for (j = 0; j <= grid_points[1]-1; j++) {
      for (i = 0; i <= grid_points[0]-1; i++) {
        for (m = 0; m < 5; m++) {
          u[k][j][i][m] = 1.0;
        }
      }
    }
  }
#endif

  //---------------------------------------------------------------------
  // first store the "interpolated" values everywhere on the grid
  //---------------------------------------------------------------------
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_START(j, cit2, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
    /*for (j = 0; j <= grid_points[1]-1; j++) {*/
      eta = (double)(j) * dnym1;
      FOR_START(i, cit3, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
      /*for (i = 0; i <= grid_points[0]-1; i++) {*/
        xi = (double)(i) * dnxm1;

        FOR_START(ix, cit4, 0, 1+1, 1, cit_step_add, RND) {
        /*for (ix = 0; ix < 2; ix++) {*/
          exact_solution((double)ix, eta, zeta, &Pface[ix][0][0]);
        }
        FOR_END(cit4);

        FOR_START(iy, cit4, 0, 1+1, 1, cit_step_add, RND) {
        /*for (iy = 0; iy < 2; iy++) {*/
          exact_solution(xi, (double)iy , zeta, &Pface[iy][1][0]);
        }
        FOR_END(cit4);

        FOR_START(iz, cit4, 0, 1+1, 1, cit_step_add, RND) {
        /*for (iz = 0; iz < 2; iz++) {*/
          exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
        }
        FOR_END(cit4);

        FOR_START(m, cit4, 0, 4+1, 1, cit_step_add, RND) {
        /*for (m = 0; m < 5; m++) {*/
          Pxi   = xi   * Pface[1][0][m] + (1.0-xi)   * Pface[0][0][m];
          Peta  = eta  * Pface[1][1][m] + (1.0-eta)  * Pface[0][1][m];
          Pzeta = zeta * Pface[1][2][m] + (1.0-zeta) * Pface[0][2][m];

          u[k][j][i][m] = Pxi + Peta + Pzeta -
                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
                          Pxi*Peta*Pzeta;
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      for (i = 0; i <= grid_points[0]-1; i++) {
        xi = (double)(i) * dnxm1;

        for (ix = 0; ix < 2; ix++) {
          exact_solution((double)ix, eta, zeta, &Pface[ix][0][0]);
        }

        for (iy = 0; iy < 2; iy++) {
          exact_solution(xi, (double)iy , zeta, &Pface[iy][1][0]);
        }

        for (iz = 0; iz < 2; iz++) {
          exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
        }

        for (m = 0; m < 5; m++) {
          Pxi   = xi   * Pface[1][0][m] + (1.0-xi)   * Pface[0][0][m];
          Peta  = eta  * Pface[1][1][m] + (1.0-eta)  * Pface[0][1][m];
          Pzeta = zeta * Pface[1][2][m] + (1.0-zeta) * Pface[0][2][m];

          u[k][j][i][m] = Pxi + Peta + Pzeta -
                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta +
                          Pxi*Peta*Pzeta;
        }
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // now store the exact values on the boundaries
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // west face
  //---------------------------------------------------------------------
  i = 0;
  xi = 0.0;
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_START(j, cit2, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
    /*for (j = 0; j <= grid_points[1]-1; j++) {*/
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // east face
  //---------------------------------------------------------------------
  i = grid_points[0]-1;
  xi = 1.0;
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_START(j, cit2, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
    /*for (j = 0; j <= grid_points[1]-1; j++) {*/
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (j = 0; j <= grid_points[1]-1; j++) {
      eta = (double)(j) * dnym1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // south face
  //---------------------------------------------------------------------
  j = 0;
  eta = 0.0;
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_START(i, cit2, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
    /*for (i = 0; i <= grid_points[0]-1; i++) {*/
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // north face
  //---------------------------------------------------------------------
  j = grid_points[1]-1;
  eta = 1.0;
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-1+1, 1, cit_step_add, RND) {
  /*for (k = 0; k <= grid_points[2]-1; k++) {*/
    zeta = (double)(k) * dnzm1;
    FOR_START(i, cit2, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
    /*for (i = 0; i <= grid_points[0]-1; i++) {*/
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 0; k <= grid_points[2]-1; k++) {
    zeta = (double)(k) * dnzm1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // bottom face
  //---------------------------------------------------------------------
  k = 0;
  zeta = 0.0;
#ifdef USE_CITERATOR
  FOR_START(j, cit1, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
  /*for (j = 0; j <= grid_points[1]-1; j++) {*/
    eta = (double)(j) * dnym1;
    FOR_START(i, cit2, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
    /*for (i =0; i <= grid_points[0]-1; i++) {*/
      xi = (double)(i) *dnxm1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (j = 0; j <= grid_points[1]-1; j++) {
    eta = (double)(j) * dnym1;
    for (i =0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) *dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // top face
  //---------------------------------------------------------------------
  k = grid_points[2]-1;
  zeta = 1.0;
#ifdef USE_CITERATOR
  FOR_START(j, cit1, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
  /*for (j = 0; j <= grid_points[1]-1; j++) {*/
    eta = (double)(j) * dnym1;
    FOR_START(i, cit2, 0, grid_points[0]-1+1, 1, cit_step_add, RND) {
    /*for (i = 0; i <= grid_points[0]-1; i++) {*/
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      FOR_START(m, cit3, 0, 4+1, 1, cit_step_add, RND) {
      /*for (m = 0; m < 5; m++) {*/
        u[k][j][i][m] = temp[m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (j = 0; j <= grid_points[1]-1; j++) {
    eta = (double)(j) * dnym1;
    for (i = 0; i <= grid_points[0]-1; i++) {
      xi = (double)(i) * dnxm1;
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
#endif // USE_CITERATOR
}


void lhsinit(double lhs[][3][5][5], int size)
{
  int i, m, n;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2;
#endif // USE_CITERATOR

  i = size;
  //---------------------------------------------------------------------
  // zero the whole left hand side for starters
  //---------------------------------------------------------------------
#ifdef USE_CITERATOR
  FOR_START(n, cit1, 0, 4+1, 1, cit_step_add, RND) {
  /*for (n = 0; n < 5; n++) {*/
    FOR_START(m, cit2, 0, 4+1, 1, cit_step_add, RND) {
    /*for (m = 0; m < 5; m++) {*/
      lhs[0][0][n][m] = 0.0;
      lhs[0][1][n][m] = 0.0;
      lhs[0][2][n][m] = 0.0;
      lhs[i][0][n][m] = 0.0;
      lhs[i][1][n][m] = 0.0;
      lhs[i][2][n][m] = 0.0;
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[0][0][n][m] = 0.0;
      lhs[0][1][n][m] = 0.0;
      lhs[0][2][n][m] = 0.0;
      lhs[i][0][n][m] = 0.0;
      lhs[i][1][n][m] = 0.0;
      lhs[i][2][n][m] = 0.0;
    }
  }
#endif // USE_CITERATOR

  //---------------------------------------------------------------------
  // next, set all diagonal values to 1. This is overkill, but convenient
  //---------------------------------------------------------------------
#ifdef USE_CITERATOR
  FOR_START(m, cit1, 0, 4+1, 1, cit_step_add, RND) {
  /*for (m = 0; m < 5; m++) {*/
    lhs[0][1][m][m] = 1.0;
    lhs[i][1][m][m] = 1.0;
  }
  FOR_END(cit1);
#else
  for (m = 0; m < 5; m++) {
    lhs[0][1][m][m] = 1.0;
    lhs[i][1][m][m] = 1.0;
  }
#endif // USE_CITERATOR
}
