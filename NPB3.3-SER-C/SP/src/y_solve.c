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
// this function performs the solution of the approximate factorization
// step in the y-direction for all five matrix components
// simultaneously. The Thomas algorithm is employed to solve the
// systems for the y-lines. Boundary conditions are non-periodic
//---------------------------------------------------------------------
void y_solve()
{
  int i, j, k, j1, j2, m;
  double ru1, fac1, fac2;
#ifdef USE_CITERATOR
  struct cit_data *cit1, *cit2, *cit3, *cit4;
#endif // USE_CITERATOR

  if (timeron) timer_start(t_ysolve);
#ifdef USE_CITERATOR
  FOR_START(k, cit1, 0, grid_points[2]-2+1, 1, cit_step_add, RND) {
  /*for (k = 1; k <= grid_points[2]-2; k++) {*/
    lhsinitj(ny2+1, nx2);

    //---------------------------------------------------------------------
    // Computes the left hand side for the three y-factors
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue
    //---------------------------------------------------------------------
    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      FOR_START(j, cit3, 0, grid_points[1]-1+1, 1, cit_step_add, RND) {
      /*for (j = 0; j <= grid_points[1]-1; j++) {*/
        ru1 = c3c4*rho_i[k][j][i];
        cv[j] = vs[k][j][i];
        rhoq[j] = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
      }
      FOR_END(cit3);

      FOR_START(j, cit3, 1, grid_points[1]-2+1, 1, cit_step_add, RND) {
      /*for (j = 1; j <= grid_points[1]-2; j++) {*/
        lhs[j][i][0] =  0.0;
        lhs[j][i][1] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1];
        lhs[j][i][2] =  1.0 + c2dtty1 * rhoq[j];
        lhs[j][i][3] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1];
        lhs[j][i][4] =  0.0;
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // add fourth order dissipation
    //---------------------------------------------------------------------
    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      j = 1;
      lhs[j][i][2] = lhs[j][i][2] + comz5;
      lhs[j][i][3] = lhs[j][i][3] - comz4;
      lhs[j][i][4] = lhs[j][i][4] + comz1;

      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz6;
      lhs[j+1][i][3] = lhs[j+1][i][3] - comz4;
      lhs[j+1][i][4] = lhs[j+1][i][4] + comz1;
    }
    FOR_END(cit2);

    FOR_START(j, cit2, 3, grid_points[1]-4+1, 1, cit_step_add, RND) {
    /*for (j = 3; j <= grid_points[1]-4; j++) {*/
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        lhs[j][i][0] = lhs[j][i][0] + comz1;
        lhs[j][i][1] = lhs[j][i][1] - comz4;
        lhs[j][i][2] = lhs[j][i][2] + comz6;
        lhs[j][i][3] = lhs[j][i][3] - comz4;
        lhs[j][i][4] = lhs[j][i][4] + comz1;
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);

    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      j = grid_points[1]-3;
      lhs[j][i][0] = lhs[j][i][0] + comz1;
      lhs[j][i][1] = lhs[j][i][1] - comz4;
      lhs[j][i][2] = lhs[j][i][2] + comz6;
      lhs[j][i][3] = lhs[j][i][3] - comz4;

      lhs[j+1][i][0] = lhs[j+1][i][0] + comz1;
      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz5;
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // subsequently, for (the other two factors
    //---------------------------------------------------------------------
    FOR_START(j, cit2, 1, grid_points[1]-2+1, 1, cit_step_add, RND) {
    /*for (j = 1; j <= grid_points[1]-2; j++) {*/
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        lhsp[j][i][0] = lhs[j][i][0];
        lhsp[j][i][1] = lhs[j][i][1] - dtty2 * speed[k][j-1][i];
        lhsp[j][i][2] = lhs[j][i][2];
        lhsp[j][i][3] = lhs[j][i][3] + dtty2 * speed[k][j+1][i];
        lhsp[j][i][4] = lhs[j][i][4];
        lhsm[j][i][0] = lhs[j][i][0];
        lhsm[j][i][1] = lhs[j][i][1] + dtty2 * speed[k][j-1][i];
        lhsm[j][i][2] = lhs[j][i][2];
        lhsm[j][i][3] = lhs[j][i][3] - dtty2 * speed[k][j+1][i];
        lhsm[j][i][4] = lhs[j][i][4];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);


    //---------------------------------------------------------------------
    // FORWARD ELIMINATION
    //---------------------------------------------------------------------
    FOR_START(j, cit2, 0, grid_points[1]-3+1, 1, cit_step_add, FWD) {
    /*for (j = 0; j <= grid_points[1]-3; j++) {*/
      j1 = j + 1;
      j2 = j + 2;
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        fac1 = 1.0/lhs[j][i][2];
        lhs[j][i][3] = fac1*lhs[j][i][3];
        lhs[j][i][4] = fac1*lhs[j][i][4];
        FOR_START(m, cit4, 0, 3, 1, cit_step_add, RND) {
        /*for (m = 0; m < 3; m++) {*/
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        FOR_END(cit4);
        lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
        lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
        FOR_START(m, cit4, 0, 3, 1, cit_step_add, RND) {
        /*for (m = 0; m < 3; m++) {*/
          rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
        }
        FOR_END(cit4);
        lhs[j2][i][1] = lhs[j2][i][1] - lhs[j2][i][0]*lhs[j][i][3];
        lhs[j2][i][2] = lhs[j2][i][2] - lhs[j2][i][0]*lhs[j][i][4];
        FOR_START(m, cit4, 0, 3, 1, cit_step_add, RND) {
        /*for (m = 0; m < 3; m++) {*/
          rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs[j2][i][0]*rhs[k][j][i][m];
        }
        FOR_END(cit4);
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      fac1 = 1.0/lhs[j][i][2];
      lhs[j][i][3] = fac1*lhs[j][i][3];
      lhs[j][i][4] = fac1*lhs[j][i][4];
      FOR_START(m, cit3, 0, 3, 1, cit_step_add, RND) {
      /*for (m = 0; m < 3; m++) {*/
        rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
      }
      FOR_END(cit3);
      lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
      lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
      FOR_START(m, cit3, 0, 3, 1, cit_step_add, RND) {
      /*for (m = 0; m < 3; m++) {*/
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
      }
      FOR_END(cit3);
      //---------------------------------------------------------------------
      // scale the last row immediately
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[j1][i][2];
      FOR_START(m, cit3, 0, 3, 1, cit_step_add, RND) {
      /*for (m = 0; m < 3; m++) {*/
        rhs[k][j1][i][m] = fac2*rhs[k][j1][i][m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors
    //---------------------------------------------------------------------
    FOR_START(j, cit2, 0, grid_points[1]-3+1, 1, cit_step_add, FWD) {
    /*for (j = 0; j <= grid_points[1]-3; j++) {*/
      j1 = j + 1;
      j2 = j + 2;
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        m = 3;
        fac1 = 1.0/lhsp[j][i][2];
        lhsp[j][i][3]    = fac1*lhsp[j][i][3];
        lhsp[j][i][4]    = fac1*lhsp[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
        lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];
        lhsp[j2][i][1]   = lhsp[j2][i][1] - lhsp[j2][i][0]*lhsp[j][i][3];
        lhsp[j2][i][2]   = lhsp[j2][i][2] - lhsp[j2][i][0]*lhsp[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsp[j2][i][0]*rhs[k][j][i][m];

        m = 4;
        fac1 = 1.0/lhsm[j][i][2];
        lhsm[j][i][3]    = fac1*lhsm[j][i][3];
        lhsm[j][i][4]    = fac1*lhsm[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
        lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];
        lhsm[j2][i][1]   = lhsm[j2][i][1] - lhsm[j2][i][0]*lhsm[j][i][3];
        lhsm[j2][i][2]   = lhsm[j2][i][2] - lhsm[j2][i][0]*lhsm[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsm[j2][i][0]*rhs[k][j][i][m];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      m = 3;
      fac1 = 1.0/lhsp[j][i][2];
      lhsp[j][i][3]    = fac1*lhsp[j][i][3];
      lhsp[j][i][4]    = fac1*lhsp[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
      lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];

      m = 4;
      fac1 = 1.0/lhsm[j][i][2];
      lhsm[j][i][3]    = fac1*lhsm[j][i][3];
      lhsm[j][i][4]    = fac1*lhsm[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
      lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];

      //---------------------------------------------------------------------
      // Scale the last row immediately
      //---------------------------------------------------------------------
      rhs[k][j1][i][3]   = rhs[k][j1][i][3]/lhsp[j1][i][2];
      rhs[k][j1][i][4]   = rhs[k][j1][i][4]/lhsm[j1][i][2];
    }
    FOR_END(cit2);


    //---------------------------------------------------------------------
    // BACKSUBSTITUTION
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    FOR_START(i, cit2, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
    /*for (i = 1; i <= grid_points[0]-2; i++) {*/
      FOR_START(m, cit3, 0, 3, 1, cit_step_add, RND) {
      /*for (m = 0; m < 3; m++) {*/
        rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3]*rhs[k][j1][i][m];
      }
      FOR_END(cit3);

      rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3]*rhs[k][j1][i][3];
      rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3]*rhs[k][j1][i][4];
    }
    FOR_END(cit2);

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    FOR_START(j, cit2, grid_points[1]-3, -1, -1, cit_step_add, FWD) {
    /*for (j = grid_points[1]-3; j >= 0; j--) {*/
      j1 = j + 1;
      j2 = j + 2;
      FOR_START(i, cit3, 1, grid_points[0]-2+1, 1, cit_step_add, RND) {
      /*for (i = 1; i <= grid_points[0]-2; i++) {*/
        FOR_START(m, cit4, 0, 3, 1, cit_step_add, RND) {
        /*for (m = 0; m < 3; m++) {*/
          rhs[k][j][i][m] = rhs[k][j][i][m] -
                            lhs[j][i][3]*rhs[k][j1][i][m] -
                            lhs[j][i][4]*rhs[k][j2][i][m];
        }
        FOR_END(cit4);

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] -
                          lhsp[j][i][3]*rhs[k][j1][i][3] -
                          lhsp[j][i][4]*rhs[k][j2][i][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] -
                          lhsm[j][i][3]*rhs[k][j1][i][4] -
                          lhsm[j][i][4]*rhs[k][j2][i][4];
      }
      FOR_END(cit3);
    }
    FOR_END(cit2);
  }
  FOR_END(cit1);
#else
  for (k = 1; k <= grid_points[2]-2; k++) {
    lhsinitj(ny2+1, nx2);

    //---------------------------------------------------------------------
    // Computes the left hand side for the three y-factors
    //---------------------------------------------------------------------

    //---------------------------------------------------------------------
    // first fill the lhs for the u-eigenvalue
    //---------------------------------------------------------------------
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (j = 0; j <= grid_points[1]-1; j++) {
        ru1 = c3c4*rho_i[k][j][i];
        cv[j] = vs[k][j][i];
        rhoq[j] = max(max(dy3+con43*ru1, dy5+c1c5*ru1), max(dymax+ru1, dy1));
      }

      for (j = 1; j <= grid_points[1]-2; j++) {
        lhs[j][i][0] =  0.0;
        lhs[j][i][1] = -dtty2 * cv[j-1] - dtty1 * rhoq[j-1];
        lhs[j][i][2] =  1.0 + c2dtty1 * rhoq[j];
        lhs[j][i][3] =  dtty2 * cv[j+1] - dtty1 * rhoq[j+1];
        lhs[j][i][4] =  0.0;
      }
    }

    //---------------------------------------------------------------------
    // add fourth order dissipation
    //---------------------------------------------------------------------
    for (i = 1; i <= grid_points[0]-2; i++) {
      j = 1;
      lhs[j][i][2] = lhs[j][i][2] + comz5;
      lhs[j][i][3] = lhs[j][i][3] - comz4;
      lhs[j][i][4] = lhs[j][i][4] + comz1;

      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz6;
      lhs[j+1][i][3] = lhs[j+1][i][3] - comz4;
      lhs[j+1][i][4] = lhs[j+1][i][4] + comz1;
    }

    for (j = 3; j <= grid_points[1]-4; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        lhs[j][i][0] = lhs[j][i][0] + comz1;
        lhs[j][i][1] = lhs[j][i][1] - comz4;
        lhs[j][i][2] = lhs[j][i][2] + comz6;
        lhs[j][i][3] = lhs[j][i][3] - comz4;
        lhs[j][i][4] = lhs[j][i][4] + comz1;
      }
    }

    for (i = 1; i <= grid_points[0]-2; i++) {
      j = grid_points[1]-3;
      lhs[j][i][0] = lhs[j][i][0] + comz1;
      lhs[j][i][1] = lhs[j][i][1] - comz4;
      lhs[j][i][2] = lhs[j][i][2] + comz6;
      lhs[j][i][3] = lhs[j][i][3] - comz4;

      lhs[j+1][i][0] = lhs[j+1][i][0] + comz1;
      lhs[j+1][i][1] = lhs[j+1][i][1] - comz4;
      lhs[j+1][i][2] = lhs[j+1][i][2] + comz5;
    }

    //---------------------------------------------------------------------
    // subsequently, for (the other two factors
    //---------------------------------------------------------------------
    for (j = 1; j <= grid_points[1]-2; j++) {
      for (i = 1; i <= grid_points[0]-2; i++) {
        lhsp[j][i][0] = lhs[j][i][0];
        lhsp[j][i][1] = lhs[j][i][1] - dtty2 * speed[k][j-1][i];
        lhsp[j][i][2] = lhs[j][i][2];
        lhsp[j][i][3] = lhs[j][i][3] + dtty2 * speed[k][j+1][i];
        lhsp[j][i][4] = lhs[j][i][4];
        lhsm[j][i][0] = lhs[j][i][0];
        lhsm[j][i][1] = lhs[j][i][1] + dtty2 * speed[k][j-1][i];
        lhsm[j][i][2] = lhs[j][i][2];
        lhsm[j][i][3] = lhs[j][i][3] - dtty2 * speed[k][j+1][i];
        lhsm[j][i][4] = lhs[j][i][4];
      }
    }


    //---------------------------------------------------------------------
    // FORWARD ELIMINATION
    //---------------------------------------------------------------------
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        fac1 = 1.0/lhs[j][i][2];
        lhs[j][i][3] = fac1*lhs[j][i][3];
        lhs[j][i][4] = fac1*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
        }
        lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
        lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
        }
        lhs[j2][i][1] = lhs[j2][i][1] - lhs[j2][i][0]*lhs[j][i][3];
        lhs[j2][i][2] = lhs[j2][i][2] - lhs[j2][i][0]*lhs[j][i][4];
        for (m = 0; m < 3; m++) {
          rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhs[j2][i][0]*rhs[k][j][i][m];
        }
      }
    }

    //---------------------------------------------------------------------
    // The last two rows in this grid block are a bit different,
    // since they for (not have two more rows available for the
    // elimination of off-diagonal entries
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      fac1 = 1.0/lhs[j][i][2];
      lhs[j][i][3] = fac1*lhs[j][i][3];
      lhs[j][i][4] = fac1*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = fac1*rhs[k][j][i][m];
      }
      lhs[j1][i][2] = lhs[j1][i][2] - lhs[j1][i][1]*lhs[j][i][3];
      lhs[j1][i][3] = lhs[j1][i][3] - lhs[j1][i][1]*lhs[j][i][4];
      for (m = 0; m < 3; m++) {
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhs[j1][i][1]*rhs[k][j][i][m];
      }
      //---------------------------------------------------------------------
      // scale the last row immediately
      //---------------------------------------------------------------------
      fac2 = 1.0/lhs[j1][i][2];
      for (m = 0; m < 3; m++) {
        rhs[k][j1][i][m] = fac2*rhs[k][j1][i][m];
      }
    }

    //---------------------------------------------------------------------
    // for (the u+c and the u-c factors
    //---------------------------------------------------------------------
    for (j = 0; j <= grid_points[1]-3; j++) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        m = 3;
        fac1 = 1.0/lhsp[j][i][2];
        lhsp[j][i][3]    = fac1*lhsp[j][i][3];
        lhsp[j][i][4]    = fac1*lhsp[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
        lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];
        lhsp[j2][i][1]   = lhsp[j2][i][1] - lhsp[j2][i][0]*lhsp[j][i][3];
        lhsp[j2][i][2]   = lhsp[j2][i][2] - lhsp[j2][i][0]*lhsp[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsp[j2][i][0]*rhs[k][j][i][m];

        m = 4;
        fac1 = 1.0/lhsm[j][i][2];
        lhsm[j][i][3]    = fac1*lhsm[j][i][3];
        lhsm[j][i][4]    = fac1*lhsm[j][i][4];
        rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
        lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
        lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
        rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];
        lhsm[j2][i][1]   = lhsm[j2][i][1] - lhsm[j2][i][0]*lhsm[j][i][3];
        lhsm[j2][i][2]   = lhsm[j2][i][2] - lhsm[j2][i][0]*lhsm[j][i][4];
        rhs[k][j2][i][m] = rhs[k][j2][i][m] - lhsm[j2][i][0]*rhs[k][j][i][m];
      }
    }

    //---------------------------------------------------------------------
    // And again the last two rows separately
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      m = 3;
      fac1 = 1.0/lhsp[j][i][2];
      lhsp[j][i][3]    = fac1*lhsp[j][i][3];
      lhsp[j][i][4]    = fac1*lhsp[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsp[j1][i][2]   = lhsp[j1][i][2] - lhsp[j1][i][1]*lhsp[j][i][3];
      lhsp[j1][i][3]   = lhsp[j1][i][3] - lhsp[j1][i][1]*lhsp[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsp[j1][i][1]*rhs[k][j][i][m];

      m = 4;
      fac1 = 1.0/lhsm[j][i][2];
      lhsm[j][i][3]    = fac1*lhsm[j][i][3];
      lhsm[j][i][4]    = fac1*lhsm[j][i][4];
      rhs[k][j][i][m]  = fac1*rhs[k][j][i][m];
      lhsm[j1][i][2]   = lhsm[j1][i][2] - lhsm[j1][i][1]*lhsm[j][i][3];
      lhsm[j1][i][3]   = lhsm[j1][i][3] - lhsm[j1][i][1]*lhsm[j][i][4];
      rhs[k][j1][i][m] = rhs[k][j1][i][m] - lhsm[j1][i][1]*rhs[k][j][i][m];

      //---------------------------------------------------------------------
      // Scale the last row immediately
      //---------------------------------------------------------------------
      rhs[k][j1][i][3]   = rhs[k][j1][i][3]/lhsp[j1][i][2];
      rhs[k][j1][i][4]   = rhs[k][j1][i][4]/lhsm[j1][i][2];
    }


    //---------------------------------------------------------------------
    // BACKSUBSTITUTION
    //---------------------------------------------------------------------
    j  = grid_points[1]-2;
    j1 = grid_points[1]-1;
    for (i = 1; i <= grid_points[0]-2; i++) {
      for (m = 0; m < 3; m++) {
        rhs[k][j][i][m] = rhs[k][j][i][m] - lhs[j][i][3]*rhs[k][j1][i][m];
      }

      rhs[k][j][i][3] = rhs[k][j][i][3] - lhsp[j][i][3]*rhs[k][j1][i][3];
      rhs[k][j][i][4] = rhs[k][j][i][4] - lhsm[j][i][3]*rhs[k][j1][i][4];
    }

    //---------------------------------------------------------------------
    // The first three factors
    //---------------------------------------------------------------------
    for (j = grid_points[1]-3; j >= 0; j--) {
      j1 = j + 1;
      j2 = j + 2;
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 3; m++) {
          rhs[k][j][i][m] = rhs[k][j][i][m] -
                            lhs[j][i][3]*rhs[k][j1][i][m] -
                            lhs[j][i][4]*rhs[k][j2][i][m];
        }

        //-------------------------------------------------------------------
        // And the remaining two
        //-------------------------------------------------------------------
        rhs[k][j][i][3] = rhs[k][j][i][3] -
                          lhsp[j][i][3]*rhs[k][j1][i][3] -
                          lhsp[j][i][4]*rhs[k][j2][i][3];
        rhs[k][j][i][4] = rhs[k][j][i][4] -
                          lhsm[j][i][3]*rhs[k][j1][i][4] -
                          lhsm[j][i][4]*rhs[k][j2][i][4];
      }
    }
  }
#endif // USE_CITERATOR
  if (timeron) timer_stop(t_ysolve);

  pinvr();
}

