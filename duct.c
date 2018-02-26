/* duct.c *************************************\
*  - Turbulent Flow in Duct modeled by LES -   *
*        two-stage time differencing           *
*  written by Andrei Chernousov Mar 16, 2001   *
*  E-mail: <andrei99@iname.com>                *
*  http://www.geocities.com/andrei_chernousov  *
\**********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <signal.h>

/* continuation flag (to be changed by SIGINT handler) */
int continFlag = 1;

/* SIGINT handler */
#pragma argsused
void termHandler(int sigType)
{
  sigType = sigType;
  if (continFlag--); /* 1 -> 0 */
  fprintf(stderr, "caught [Ctrl+C], terminating...");
}

/*--- Combinations of the ratio of specific heats ---*/
#define K        1.4
#define K_1      0.4
#define _K_1     (-0.4)
#define _05K     0.7
#define K_K_1    3.5
#define K_1_2    0.2
#define _1_2_K_1 1.25

#define Cv      717.75  /* specific constant volume heat capacity, J/(kg*K) */
#define R_VOZD  287.10  /* gas constant for air */
#define KR_VOZD 401.94  /* R_VOZD * K */

#define Fx1 xU1 /* aliases for x-fluxes */
#define Fx2 xU2
#define Fx3 xU3
#define Fx4 xU4
#define Fx5 xU5

#define Fy1 yU1 /* aliases for y-fluxes */
#define Fy2 yU2
#define Fy3 yU3
#define Fy4 yU4
#define Fy5 yU5

#define Fz1 zU1 /* aliases for z-fluxes */
#define Fz2 zU2
#define Fz3 zU3
#define Fz4 zU4
#define Fz5 zU5

/*=== Global variables ===*/
float
  ***U1_, ***U2_, ***U3_, ***U4_, ***U5_, /* substitutional pointers */
  /* Cons. variables at the cells' centroids */
  ***U1,  ***U2,  ***U3,  ***U4,  ***U5,
  ***U1p, ***U2p, ***U3p, ***U4p, ***U5p,
  /* Cons. variables at cell boundaries in x-direction */
  ***xU1, ***xU2, ***xU3, ***xU4, ***xU5,
  ***U1x, ***U2x, ***U3x, ***U4x, ***U5x,
  /* Cons. variables at cell boundaries in y-direction */
  ***yU1, ***yU2, ***yU3, ***yU4, ***yU5,
  ***U1y, ***U2y, ***U3y, ***U4y, ***U5y,
  /* Cons. variables at cell boundaries in z-direction */
  ***zU1, ***zU2, ***zU3, ***zU4, ***zU5,
  ***U1z, ***U2z, ***U3z, ***U4z, ***U5z,

  R0, /* initial density */
  R, U, V, W, P, C,       /* vector of primitive flow parameters */
  u1, u2, u3, u4, u5,     /* conservative flow variables         */
  deltaX, deltaY, deltaZ, /* cell spacings                       */
  deltaT,                 /* step in time                        */
  deltaT_X, deltaT_Y, deltaT_Z, /* convenient ratios             */
  _deltaX, _deltaY, _deltaZ,
  _4deltaX, _4deltaY, _4deltaZ,
  _2deltaX, _2deltaY, _2deltaZ;

char
  Stage, /* indicator of the current stage */
  str[80]; /* buffer string */

unsigned
  LEN, HIG, DEP, /* numbers o cells in x-, y-, z-directions, for 3d box */
  LENN, HIGG, DEPP,
  Answer,  /* flag of continuing solution */
  i,  j,  k, i__,  j__,  k__, /* indices */
  l,

  mem,     /* total amount of dynamic memory required (in bytes) */
  step,    /* current time step */
  f_step,  /* frame taking step */
  numstep; /* overall prescribed number of time steps */

FILE *pF; /* pointer to a file structure */

/*--- For monotone piecevise parabolic reconstruction ---*/
/* constant coefficients */
#define BB 4.000000000000000
#define C1 0.333333333333333
#define C2 0.166666666666667

/***************
* minmod(x, y) * - macro: limiting function of Chakravarthy, Osher
***************/
/*
	- uses register variable kk defined in Reconstruction()
	- don't use expressions as arguments of this macro!
	- beware of other side effects
*/
#define minmod(x, y) ( (x*y <= 0.) ? 0.: \
  ( (kk=(y>0.)?y*BB:-y*BB, x>0.) ? ((x<kk)?x:kk) : ((-x<kk)?x:-kk) ) )

/***************
* minmod(x, y) * - function: limiting function of Chakravarthy, Osher (slow!)
***************/
/*
float ay;
float minmod(float x, float y)
{
  if (x * y <= 0.)
    return 0.;
  ay = ((y *= BB) > 0.) ? y: -y;
  if (x > 0.)
    return (x < ay) ? x : ay;
  return (-x < ay) ? x: -ay;
} */ /* end minmod() */


/* finite differences of conservative variables */
float
  m1, m2, m3, m4, m5,
  w1, w2, w3, w4, w5,
  /* finite differences of characteristic variables */
  M1, M2, M3, M4, M5,
  W1, W2, W3, W4, W5,
  /* modified finite differences of characteristic variables */
  _M1, _M2, _M3, _M4, _M5,
  _W1, _W2, _W3, _W4, _W5,
  /* transformation matrix [S-1] */
  S11, S12, S13,
  S23,
  S41, S42,
  S51, S52,
  /* inverted transformation matrix [S-1] */
  _S23, _S24, _S25,
  _S33,
  _S43,
  _S53, _S54, _S55,
  /* subsidiary */
  aa, bb, cc, dd, ee, ff, gg, hh, ll, mm, nn,
  u_hh, v_hh, w_hh, u_ee, v_ee, w_ee,
  k1, k2, k3, k4, k5,

  /*--- For LCS Riemann solver ---*/
  P_, _P, R_, _R, U_, _U, V_, _V, W_, _W, C_, _C, _T, T_,
  C_p, C_m, C_o, C_p_, C_m_, C_o_,
  _jo, _jp, _jm, jo_, jp_, jm_,
  al_p, al_m, al_o,
  J_p, J_m, J_o,

  /*--- For evaluation of the viscous fluxes ---*/
  /* matrix of velocity derivatives */
  du_dx, dv_dx, dw_dx,
  du_dy, dv_dy, dw_dy,
  du_dz, dv_dz, dw_dz,
  /* velocity components at surrounding cells */
  iu, ui, ju, uj, ku, uk,
  iv, vi, jv, vj, kv, vk,
  iw, wi, jw, wj, kw, wk,
  iU, Ui, jU, Uj, kU, Uk,
  iV, Vi, jV, Vj, kV, Vk,
  iW, Wi, jW, Wj, kW, Wk,
  /* mean filrered velocity components and density at the interfaces */
  uu, vv, ww, rr,
  /* viscous stress tensor */
  sigma_xx, sigma_xy, sigma_xz,
  sigma_yx, sigma_yy, sigma_yz,
  sigma_zx, sigma_zy, sigma_zz,
  /* heat flux vector */
  q_x, q_y, q_z,
  /* transport coefficients */
  /* molecular */
  mu_L, lambda_L, Pr_L,
  /* subgrid-scale */
  mu_T, Pr_T, cp_Pr_T,
  /* effective */
  mu_E, lambda_E,
  /* Smagorinsky constant and a complex CsDD = Cs * delta * delta, where delta = (dX*dY*dZ)~0.33333 */
  Cs, CsDD,
  /* arrays of corrected Smagorinsky complex */
  **CsDDx, **CsDDy, **CsDDz,
  /* invariant of deforming rate tensor */
  _S_,
  /* Cartesian coordinates */
  x, y, z,
  /* distances in wall units, friction velocity, mean wall stress */
  d_plus, u_tau, tau_w,
  /* accceleration of mass force in source term */
  g, gdeltaT,
  /* array of probes */
  *probes,
  M, I, E,    /* mean fluxes */
  tauY, tauZ; /* mean wall stresses */

   /* files to output data */
FILE *pFprobest, *pFprobesb;


/************
*  ARRAY2D  *   Memory allocation procedure for 2D arrays
************/
float **Array2D(unsigned columns, unsigned rows)
{
  float **x;
  unsigned i;

  if ((x = (float **)malloc(columns*sizeof(float*))) == NULL) {
    fprintf(stderr, "can't allocate memory.\n");
    exit(-1);
  }
  for (i = 0; i < columns; i++)
    if ((x[i] = (float *)malloc(rows*sizeof(float))) == NULL) {
      fprintf(stderr, "can't allocate memory.\n");
      exit(-1);
    }

  return x;

} /* end Array2D() */


/************
*  ARRAY3D  *   Memory allocation procedure for 3D arrays
************/
float ***Array3D(unsigned columns, unsigned rows, unsigned floors)
{
  float ***x;
  unsigned i, j;

  if ((x = (float ***)malloc(columns*sizeof(float**))) == NULL) {
    fprintf(stderr, "can't allocate memory.\n");
    exit(-1);
  }
  for (i = 0; i < columns; i++) {
    if ((x[i] = (float **)malloc(rows*sizeof(float*))) == NULL) {
      fprintf(stderr, "can't allocate memory.\n");
      exit(-1);
    }
    for (j = 0; j < rows; j++)
      if ((x[i][j] = (float *)malloc(floors*sizeof(float))) == NULL) {
	fprintf(stderr, "can't allocate memory.\n");
	exit(-1);
      }
  }

  return x;

} /* end Array3D() */


/***************
*  INITIALIZE  *    Performs necessary initializations
***************/
void Initialize(void)
{
  /*--- Loading data from "duct.ini" file ---*/
  if ((pF = fopen("duct.ini", "r")) == NULL) {
    fprintf(stderr, "can't open \"duct.ini\".\n");
    exit(-1);
  }
  i = 0;
  do{
    if (fgets(str, 80, pF) != NULL) {
      switch (i) {
      case  0: printf("LEN    = %d\n", LEN    = atoi(str)); break;
      case  1: printf("HIG    = %d\n", HIG    = atoi(str)); break;
      case  2: printf("DEP    = %d\n", DEP    = atoi(str)); break;
      case  3: printf("deltaX = %f\n", deltaX = atof(str)); break;
      case  4: printf("deltaY = %f\n", deltaY = atof(str)); break;
      case  5: printf("deltaZ = %f\n", deltaZ = atof(str)); break;
      case  6: printf("deltaT = %f\n", deltaT = atof(str)); break;
      case  7: printf("Answer = %d\n", Answer = atoi(str)); break;
      case  8: printf("numstep= %d\n", numstep= atoi(str)); break;
      case  9: printf("mu_L   = %f\n", mu_L   = atof(str)); break;
      case 10: printf("Pr_L   = %f\n", Pr_L   = atof(str)); break;
      case 11: printf("Pr_T   = %f\n", Pr_T   = atof(str)); break;
      case 12: printf("Cs     = %f\n", Cs     = atof(str)); break;
      case 13: printf("g      = %f\n", g      = atof(str)); break;
      case 14: printf("tau_w  = %f\n", tau_w  = atof(str)); break;
      case 15: printf("R0     = %f\n", R0     = atof(str)); break;
      case 16: printf("f_step = %d\n", f_step = atoi(str)); break;
      default: break;
      } /* end switch */
      i++;
    } /* end if() */
    else {
      fprintf(stderr, "error reading \"duct.ini\".\n");
      exit(-1);
    }
  } while (i < 17);

  fclose( pF );

  /*--- other initialisatons ---*/
  /* complexes with deltas */
  deltaT_X = deltaT / deltaX;     deltaT_Y = deltaT / deltaY;     deltaT_Z = deltaT / deltaZ;
  _deltaX = 1. / deltaX;          _deltaY  = 1. / deltaY;          _deltaZ = 1. / deltaZ;
  _4deltaX = 1. / ( 4 * deltaX ); _4deltaY = 1. / ( 4 * deltaY ); _4deltaZ = 1. / ( 4 * deltaZ );
  _2deltaX = _4deltaX + _4deltaX; _2deltaY = _4deltaY + _4deltaY; _2deltaZ = _4deltaZ + _4deltaZ;
  gdeltaT = g * deltaT;
  /* box sizes + 1 */
  LENN = LEN + 1; HIGG = HIG + 1; DEPP = DEP + 1;
  /* constant-value lambda_L */
  lambda_L = mu_L * Cv * K / Pr_L;
  /* complex cp_Pr_T */
  cp_Pr_T = Cv * K / Pr_T;
  /* complex Cd * delta * delta in Smagorinsky model */
  CsDD = Cs * pow(deltaX * deltaY * deltaZ, 0.6666666666667);
  /* friction velocity */
  u_tau = sqrt(tau_w / R0);
  /* x-direction near-wall correction */
  CsDDx = Array2D(HIG, DEP);
  for (j = 0; j < HIG; j++) {
    for (k = 0; k < DEP; k++) {
      y = (j < HIG/2) ? (j+0.5) * deltaY : (HIG-j-0.5) * deltaY;
      z = (k < DEP/2) ? (k+0.5) * deltaZ : (DEP-k-0.5) * deltaZ;
      dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
      d_plus = R0 * u_tau * dd / mu_L;
      CsDDx[j][k] = CsDD * (1. - exp( - (d_plus/25) * (d_plus/25) * (d_plus/25)));
    }
  }
  /* y-direction near-wall correction */
  CsDDy = Array2D(HIGG, DEP);
  for (j = 0; j < HIGG; j++) {
    for (k = 0; k < DEP; k++) {
      y = (j <= HIG/2) ? j       * deltaY : (HIG-j)     * deltaY;
      z = (k <  DEP/2) ? (k+0.5) * deltaZ : (DEP-k-0.5) * deltaZ;
      dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
      d_plus = R0 * u_tau * dd / mu_L;
      CsDDy[j][k] = CsDD * (1. - exp( - (d_plus/25) * (d_plus/25) * (d_plus/25)));
    }
  }
  /* z-direction near-wall correction */
  CsDDz = Array2D( HIG, DEPP );
  for( j = 0; j < HIG; j++ ) {
    for( k = 0; k < DEPP; k++ ) {
      y = (j  < HIG/2) ? (j+0.5) * deltaY : (HIG-j-0.5) * deltaY;
      z = (k <= DEP/2) ? k       * deltaZ : (DEP-k)     * deltaZ;
      dd = 2. * y * z / (y + z + sqrt(y * y + z * z));
      d_plus = R0 * u_tau * dd / mu_L;
      CsDDz[j][k] = CsDD * (1. - exp( - (d_plus/25) * (d_plus/25) * (d_plus/25)));
    }
  }

  /*--- Allocation of dybamic memory ---*/
  {
    /* temporal array of pointers */
    float ***fp[10];
    /* total amount of memory required, for gas dynamics 5 */
    mem
      = 10 * ( (unsigned long)LEN + 2 ) * ( (unsigned long)HIG + 2 ) * ( (unsigned long)DEP + 2 )
      + 10 * ( (unsigned long)LEN + 1 ) *   (unsigned long)HIG	     *   (unsigned long)DEP
      + 10 *   (unsigned long)LEN       * ( (unsigned long)HIG + 1 ) *   (unsigned long)DEP
      + 10 *   (unsigned long)LEN	*   (unsigned long)HIG	     * ( (unsigned long)DEP + 1 );
    printf( "%d bytes of dynamic memory required...", mem*sizeof(float) );
    /* for quantities in cells */
    for (i = 0; i < 10; i++)
      fp[i] = Array3D(LEN+2, HIG+2, DEP+2);
    U1 = fp[0], U2 = fp[1], U3 = fp[2], U4 = fp[3], U5 = fp[4];
    U1p= fp[5], U2p= fp[6], U3p= fp[7], U4p= fp[8], U5p= fp[9];
    /* for x-fluxes */
    for (i = 0; i < 10; i++)
      fp[i] = Array3D(LENN, HIG, DEP);
    xU1 = fp[0], xU2 = fp[1], xU3 = fp[2], xU4 = fp[3], xU5 = fp[4];
    U1x = fp[5], U2x = fp[6], U3x = fp[7], U4x = fp[8], U5x = fp[9];
    /* for y-fluxes */
    for (i = 0; i < 10; i++)
      fp[i] = Array3D(LEN, HIGG, DEP);
    yU1 = fp[0], yU2 = fp[1], yU3 = fp[2], yU4 = fp[3], yU5 = fp[4];
    U1y = fp[5], U2y = fp[6], U3y = fp[7], U4y = fp[8], U5y = fp[9];
    /* for z-fluxes */
    for (i = 0; i < 10; i++)
      fp[i] = Array3D(LEN, HIG, DEPP);
    zU1 = fp[0], zU2 = fp[1], zU3 = fp[2], zU4 = fp[3], zU5 = fp[4];
    U1z = fp[5], U2z = fp[6], U3z = fp[7], U4z = fp[8], U5z = fp[9];
    printf(" allocated!\n");
  }

  /* array of probes */
  if ((probes = (float*)malloc(sizeof(float)*HIG)) == NULL) {
    fprintf(stderr, "can't allocate memory for probes.\n");
    exit(-1);
  }

  /*--- Initial conditions ---*/
  if (Answer) { /* restore solution from backup */
    if ((pF = fopen("backup.dat", "rb")) == NULL) {
      fprintf(stderr, "can't open \"backup.dat\" to read.\n");
      exit(-1);
    }
    for (i = 1; i < LENN; i++) {
      for (j = 1; j < HIGG; j++) {
	fread(&U1[i][j][1], DEP * sizeof(U1[i][j][1]), 1, pF);
	fread(&U2[i][j][1], DEP * sizeof(U2[i][j][1]), 1, pF);
	fread(&U3[i][j][1], DEP * sizeof(U3[i][j][1]), 1, pF);
	fread(&U4[i][j][1], DEP * sizeof(U4[i][j][1]), 1, pF);
	fread(&U5[i][j][1], DEP * sizeof(U5[i][j][1]), 1, pF);
      }
    }
    fread(&step, sizeof(step), 1, pF);
    fclose(pF);
  }
  else {
    /* start new simulation */
    for (i = 1; i < LENN; i++) {
      U = 187.08; V = W = 0.;
      P = 100000., R = R0;
      for (j = 1; j < HIGG; j++) {
	for (k = 1; k < DEPP; k++) {
	  if (j > HIG/5 && j <= 4*HIG/5 && k > DEP/5 && k <= 4*DEP/5) {
	    if (i <= LEN/4)          { V =  65.; W =  50.; }
	    else if (i <= 2 * LEN/4) { V =  40.; W = -50.; }
	    else if (i <= 3 * LEN/4) { V = -65.; W =  35.; }
	    else                     { V = -40.; W = -35.; }
	  }
	  else
	    V = W = 0.;
	  U1[i][j][k] = R;
	  U2[i][j][k] = R * U;
	  U3[i][j][k] = R * V;
	  U4[i][j][k] = R * W;
	  U5[i][j][k] = P / K_1 + 0.5 * R * (U * U + V * V + W * W);
	}
      }
    }
    step = 0;  /* zero init step number at the beginning */
  } /* end else */

  Stage = 1;
  U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;

  /*--- open "probes.txt" file in text mode ---*/
  if ((pFprobest = fopen("probes.txt", "a")) == NULL) {
    fprintf(stderr, "can't open \"probes.txt\".\n");
    exit(-1);
  }

  /*--- open "probes.bin" file in binary mode ---*/
  if ((pFprobesb = fopen("probes.bin", "ab")) == NULL) {
    fprintf(stderr, "can't open \"probes.bin\".\n");
    exit(-1);
  }

} /* end Initialize( ) */


/***********************
* BOUNCONDINGHOSTCELLS *   BC in the cells out of the domain
***********************/
void BounCondInGhostCells( void )
{
  /*register*/ unsigned i, j, k;

  /* x */
  for( j = 1; j < HIGG; j++ ) {
    for( k = 1; k < DEPP; k++ ) {
      /* left side - PERIODIC INFLOW */
      U1_[0][j][k] = U1_[LEN][j][k];
      U2_[0][j][k] = U2_[LEN][j][k];
      U3_[0][j][k] = U3_[LEN][j][k];
      U4_[0][j][k] = U4_[LEN][j][k];
      U5_[0][j][k] = U5_[LEN][j][k];
      /* right side - PERIODIC OUTFLOW */
      U1_[LENN][j][k] = U1_[1][j][k];
      U2_[LENN][j][k] = U2_[1][j][k];
      U3_[LENN][j][k] = U3_[1][j][k];
      U4_[LENN][j][k] = U4_[1][j][k];
      U5_[LENN][j][k] = U5_[1][j][k];
    } /* end for */
  } /* end for */
  /* y */
  for( i = 0; i <= LENN; i++ ) {
    for( k = 1; k < DEPP; k++ ) {
      /* bottom side - NO SLIP */
      U1_[i][0][k] =   U1_[i][1][k];
      U2_[i][0][k] = - U2_[i][1][k];
      U3_[i][0][k] = - U3_[i][1][k];
      U4_[i][0][k] = - U4_[i][1][k];
      U5_[i][0][k] =   U5_[i][1][k];
      /* top side - NO SLIP */
      U1_[i][HIGG][k] =   U1_[i][HIG][k];
      U2_[i][HIGG][k] = - U2_[i][HIG][k];
      U3_[i][HIGG][k] = - U3_[i][HIG][k];
      U4_[i][HIGG][k] = - U4_[i][HIG][k];
      U5_[i][HIGG][k] =   U5_[i][HIG][k];
    }
  }
  /* z */
  for( i = 0; i <= LENN; i++ ) {
    for( j = 1; j < HIGG; j++ ) {
      /* back side - NO SLIP */
      U1_[i][j][0] =   U1_[i][j][1];
      U2_[i][j][0] = - U2_[i][j][1];
      U3_[i][j][0] = - U3_[i][j][1];
      U4_[i][j][0] = - U4_[i][j][1];
      U5_[i][j][0] =   U5_[i][j][1];
      /* front side - NO SLIP */
      U1_[i][j][DEPP] =   U1_[i][j][DEP];
      U2_[i][j][DEPP] = - U2_[i][j][DEP];
      U3_[i][j][DEPP] = - U3_[i][j][DEP];
      U4_[i][j][DEPP] = - U4_[i][j][DEP];
      U5_[i][j][DEPP] =   U5_[i][j][DEP];
    }
  }
  /* along x */
  for( i = 0; i <= LENN; i++ ) { /* NO SLIP */
    U1_[i][0][0] =   U1_[i][1][1];
    U2_[i][0][0] = - U2_[i][1][1];
    U3_[i][0][0] = - U3_[i][1][1];
    U4_[i][0][0] = - U4_[i][1][1];
    U1_[i][0][DEPP] =   U1_[i][1][DEP];
    U2_[i][0][DEPP] = - U2_[i][1][DEP];
    U3_[i][0][DEPP] = - U3_[i][1][DEP];
    U4_[i][0][DEPP] = - U4_[i][1][DEP];
    U1_[i][HIGG][0] =   U1_[i][HIG][1];
    U2_[i][HIGG][0] = - U2_[i][HIG][1];
    U3_[i][HIGG][0] = - U3_[i][HIG][1];
    U4_[i][HIGG][0] = - U4_[i][HIG][1];
    U1_[i][HIGG][DEPP] =   U1_[i][HIG][DEP];
    U2_[i][HIGG][DEPP] = - U2_[i][HIG][DEP];
    U3_[i][HIGG][DEPP] = - U3_[i][HIG][DEP];
    U4_[i][HIGG][DEPP] = - U4_[i][HIG][DEP];
  }

} /* end BounCondInGhostCells() */


/*******************
*  RECONSTRUCTION  *   Piecewise-parabolic reconstruction (higher order in space)
*******************/
void Reconstruction(void)
{
/*register*/ float kk;
/*register*/ unsigned i, j, k, _i, _j, _k, i_, j_, k_;

  for( i = 1, _i = 0, i_ = 2; i < LENN; i++, _i++, i_++ ) {
    for( j = 1, _j = 0, j_ = 2; j < HIGG; j++, _j++, j_++ ) {
      for( k = 1, _k = 0, k_ = 2; k < DEPP; k++, _k++, k_++ ) {
	/*---  X  ---*/
	/* at cell centroid */
	u1 = U1_[i][j][k];
	u2 = U2_[i][j][k];
	u3 = U3_[i][j][k];
	u4 = U4_[i][j][k];
	u5 = U5_[i][j][k];
	R =   1. / u1;   /* 1 / R  */
	U =   u2 * R;
	V =   u3 * R;
	W =   u4 * R;
	P = ( u5 - 0.5 * ( u2 * u2 + u3 * u3 + u4 * u4 ) * R ) * K_1;
	C = sqrt( K * P * R );
				/* evaluation of matrices */
	ff = U*U+V*V+W*W;   /*                       */
	aa = K_1_2*ff;      /* (K-1)*(U*U+V*V+W*W)/2 */
	bb = _K_1*U;        /* - (K_1) * U           */
	cc = _K_1*V;        /* - (K_1) * V           */
	dd = _K_1*W;        /* - (K_1) * W           */
	ll = 1./(C+C);      /* 1/(2.*C)              */
	ee = ll/C;          /* 1./(2*C*C)            */
	gg = - ff * ee;     /* -(U*U+V*V+W*W)/(2*C*C)*/
	hh = - ( ee + ee ); /* - 1 / (C*C)           */
	mm = U*ll;          /* U/(2.*C)              */
				/* differencing of the conservative variables*/
	/* forvard */
	m1 = U1_[i_][j][k] - u1;
	m2 = U2_[i_][j][k] - u2;
	m3 = U3_[i_][j][k] - u3;
	m4 = U4_[i_][j][k] - u4;
	m5 = U5_[i_][j][k] - u5;
	/* backward */
	w1 = u1 - U1_[_i][j][k];
	w2 = u2 - U2_[_i][j][k];
	w3 = u3 - U3_[_i][j][k];
	w4 = u4 - U4_[_i][j][k];
	w5 = u5 - U5_[_i][j][k];
	/* transformation matrix [S] */
	S11 = aa - C*C;
	S41 = aa - (kk=C*U);   /* CU */
	S51 = aa + kk;
	S12 = S42 = S52 = bb;
	S42 += C;
	S52 -= C;
	/* finite differences of characteristic variables dW = S * dU */
	/* forward */
	kk = K_1 * m5;
	k2 = cc * m3;
	k3 = dd * m4;
	M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
	M2 =  -V*m1 		 +   m3            ;
	M3 =  -W*m1 		 		+ 	m4     ;
	M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
	M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
	/* backward */
	kk = K_1 * w5;
	k2 = cc * w3;
	k3 = dd * w4;
	W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
	W2 =  -V*w1 		 +   w3 	       ;
	W3 =  -W*w1 				+   w4     ;
	W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
	W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
	/* modified characteristic differences _W */
	_M1 = minmod( M1, W1 ); /* the best way */
	_M2 = minmod( M2, W2 );
	_M3 = minmod( M3, W3 );
	_M4 = minmod( M4, W4 );
	_M5 = minmod( M5, W5 );
	_W1 = minmod( W1, M1 );
	_W2 = minmod( W2, M2 );
	_W3 = minmod( W3, M3 );
	_W4 = minmod( W4, M4 );
	_W5 = minmod( W5, M5 );
	/* inverted transformation matrix [S-1] */
	/*_S21 = */ u_hh = U*hh; /*_S31 =*/ v_hh = V*hh; /*_S41 =*/ w_hh = W*hh;
	u_ee = U*ee; /*_S34 =*/ v_ee = V*ee; /*_S44 =*/ w_ee = W*ee;
	_S24 = u_ee + ll; _S25 = u_ee - ll;
	_S54 = _S55 = ( nn = - 0.5 * gg ) + _1_2_K_1;
	_S54 += mm; _S55 -= mm;
	/* parameters' vector at the rightmost interface */
	k1 = _M1 * C1 + _W1 * C2;
	k2 = _M2 * C1 + _W2 * C2;
	k3 = _M3 * C1 + _W3 * C2;
	k4 = _M4 * C1 + _W4 * C2;
	k5 = _M5 * C1 + _W5 * C2;
	kk = k4 + k5;
	xU1[ i][_j][_k] = u1 +   hh * k1 + 						 	ee  * kk;
	xU2[ i][_j][_k] = u2 + u_hh * k1 + 						   _S24 * k4 + _S25 * k5;
	xU3[ i][_j][_k] = u3 + v_hh * k1 + 		  k2 			 + v_ee * kk;
	xU4[ i][_j][_k] = u4 + w_hh * k1 + 					  k3 + w_ee * kk;
	xU5[ i][_j][_k] = u5 +   gg * k1 +    V * k2 +    W * k3 + _S54 * k4 + _S55 * k5;
	/* ... leftmost interface */
	k1 = _M1 * C2 + _W1 * C1;
	k2 = _M2 * C2 + _W2 * C1;
	k3 = _M3 * C2 + _W3 * C1;
	k4 = _M4 * C2 + _W4 * C1;
	k5 = _M5 * C2 + _W5 * C1;
	kk =  k4 + k5;
	U1x[_i][_j][_k] = u1 -   hh * k1 - 						     ee * kk;
	U2x[_i][_j][_k] = u2 - u_hh * k1 - 						   _S24 * k4 - _S25 * k5;
	U3x[_i][_j][_k] = u3 - v_hh * k1 - 		  k2 			 - v_ee * kk;
	U4x[_i][_j][_k] = u4 - w_hh * k1 - 					  k3 - w_ee * kk;
	U5x[_i][_j][_k] = u5 -   gg * k1 -    V * k2 -    W * k3 - _S54 * k4 - _S55 * k5;
	/*---  Y  ---*/
	/* evaluation of matrices */
	kk = bb; bb = cc; cc = dd; dd = kk;
	mm = V*ll;          /* V/(2.*C) */
	/* differencing of the conservative variables */
	/* forvard */
	m1 = U1_[i][j_][k] - u1;
	m2 = U3_[i][j_][k] - u3; /* U3[i][j+1][k] ( U->V, V->W, W->U : 3-4-2 ) */
	m3 = U4_[i][j_][k] - u4; /* U4[i][j+1][k] */
	m4 = U2_[i][j_][k] - u2; /* U2[i][j+1][k] */
	m5 = U5_[i][j_][k] - u5;
	/* backward */
	w1 = u1 - U1_[i][_j][k];
	w2 = u3 - U3_[i][_j][k]; /* U3[i][j-1][k] ( U->V, V->W, W->U : 3-4-2 ) */
	w3 = u4 - U4_[i][_j][k]; /* U4[i][j-1][k] */
	w4 = u2 - U2_[i][_j][k]; /* U2[i][j-1][k] */
	w5 = u5 - U5_[i][_j][k];
	/* transformation matrix [S] */
	S41 = aa - (kk=C*V);    /* CV */
	S51 = aa + kk;
	S12 = S42 = S52 = bb;
	S42 += C;
	S52 -= C;
	/* finite differences of characteristic variables dW = S * dU */
	/* forward */
	kk = K_1 * m5;
	k2 = cc * m3;
	k3 = dd * m4;
	M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
	M2 =  -W*m1 		 +   m3              ;  /* -W*m1 ! */
	M3 =  -U*m1 		 		+ 	m4          ;  /* -U*m1 ! */
	M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
	M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
	/* backward */
	kk = K_1 * w5;
	k2 = cc * w3;
	k3 = dd * w4;
	W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
	W2 =  -W*w1 		 +   w3 	       ;  /* -W*m1 ! */
	W3 =  -U*w1 				+   w4     ;  /* -U*m1 ! */
	W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
	W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
	/* modified characteristic differences _W */
	_M1 = minmod( M1, W1 );
	_M2 = minmod( M2, W2 );
	_M3 = minmod( M3, W3 );
	_M4 = minmod( M4, W4 );
	_M5 = minmod( M5, W5 );
	_W1 = minmod( W1, M1 );
	_W2 = minmod( W2, M2 );
	_W3 = minmod( W3, M3 );
	_W4 = minmod( W4, M4 );
	_W5 = minmod( W5, M5 );
	/* inverted transformation matrix [S-1] */
	/*_S21 =  v_hh;   _S31 =  w_hh; _S41 =  u_hh; */
	/* _S34 =  w_ee; _S44 =  u_ee; */
	_S24 = v_ee + ll; _S25 = v_ee - ll;
	_S54 = _S55 = nn + _1_2_K_1;
	_S54 += mm; _S55 -= mm;
	/* parameters' vector at the upper interface */
	k1 = _M1 * C1 + _W1 * C2;
	k2 = _M2 * C1 + _W2 * C2;
	k3 = _M3 * C1 + _W3 * C2;
	k4 = _M4 * C1 + _W4 * C2;
	k5 = _M5 * C1 + _W5 * C2;
	kk = k4 + k5;
	yU1[_i][ j][_k] = u1 +   hh * k1 + 						 	ee  * kk;
	yU3[_i][ j][_k] = u3 + v_hh * k1 + 						   _S24 * k4 + _S25 * k5;
	yU4[_i][ j][_k] = u4 + w_hh * k1 + 		  k2 			 + w_ee * kk;
	yU2[_i][ j][_k] = u2 + u_hh * k1 + 					  k3 + u_ee * kk;
	yU5[_i][ j][_k] = u5 +   gg * k1 +    W * k2 +    U * k3 + _S54 * k4 + _S55 * k5;
	/* ...lower interface */
	k1 = _M1 * C2 + _W1 * C1;
	k2 = _M2 * C2 + _W2 * C1;
	k3 = _M3 * C2 + _W3 * C1;
	k4 = _M4 * C2 + _W4 * C1;
	k5 = _M5 * C2 + _W5 * C1;
	kk =  k4 + k5;
	U1y[_i][_j][_k] = u1 -   hh * k1 - 						     ee * kk;
	U3y[_i][_j][_k] = u3 - v_hh * k1 - 						   _S24 * k4 - _S25 * k5;
	U4y[_i][_j][_k] = u4 - w_hh * k1 - 		  k2 			 - w_ee * kk;
	U2y[_i][_j][_k] = u2 - u_hh * k1 - 					  k3 - u_ee * kk;
	U5y[_i][_j][_k] = u5 -   gg * k1 -    W * k2 -    U * k3 - _S54 * k4 - _S55 * k5;
	/*--- Z ---*/
	/* evaluation of matrices */
	kk = bb; bb = cc; cc = dd; dd = kk;
	mm = W*ll;          /* W/(2.*C) */
	/* differencing of the conservative variables */
	/* forward */
	m1 = U1_[i][j][k_] - u1;
	m2 = U4_[i][j][k_] - u4; /* U4[i][j][k+1] ( V->W, W->U, U->V : 4-2-3 ) */
	m3 = U2_[i][j][k_] - u2; /* U2[i][j][k+1] */
	m4 = U3_[i][j][k_] - u3; /* U3[i][j][k+1] */
	m5 = U5_[i][j][k_] - u5;
	/* backward */
	w1 = u1 - U1_[i][j][_k];
	w2 = u4 - U4_[i][j][_k]; /* U4[i][j][k-1] ( V->W, W->U, U->V : 4-2-3 ) */
	w3 = u2 - U2_[i][j][_k]; /* U2[i][j][k-1] */
	w4 = u3 - U3_[i][j][_k]; /* U3[i][j][k-1] */
	w5 = u5 - U5_[i][j][_k];
	/* transformation matrix [S] */
	S41 = aa - (kk=C*W);    /* CW */
	S51 = aa + kk;
	S12 = S42 = S52 = bb;
	S42 += C;
	S52 -= C;
	/* finite differences of characteristic variables dW = S * dU */
	/* forward */
	kk = K_1 * m5;
	k2 = cc * m3;
	k3 = dd * m4;
	M1 = S11*m1 + S12*m2 +   k2 +   k3 + kk;
	M2 =  -U*m1 		 +   m3            ;  /* -U*m1 ! */
	M3 =  -V*m1 		 		+ 	m4     ;  /* -V*m1 ! */
	M4 = S41*m1 + S42*m2 +   k2 +   k3 + kk;
	M5 = S51*m1 + S52*m2 +   k2 +   k3 + kk;
	/* backward */
	kk = K_1 * w5;
	k2 = cc * w3;
	k3 = dd * w4;
	W1 = S11*w1 + S12*w2 +   k2 +   k3 + kk;
	W2 =  -U*w1 		 +   w3 	       ;   /* -U*m1 ! */
	W3 =  -V*w1 				+   w4     ;  /* -V*m1 ! */
	W4 = S41*w1 + S42*w2 +   k2 +   k3 + kk;
	W5 = S51*w1 + S52*w2 +   k2 +   k3 + kk;
	/* modified characteristic differences _W */
	_M1 = minmod( M1, W1 );
	_M2 = minmod( M2, W2 );
	_M3 = minmod( M3, W3 );
	_M4 = minmod( M4, W4 );
	_M5 = minmod( M5, W5 );
	_W1 = minmod( W1, M1 );
	_W2 = minmod( W2, M2 );
	_W3 = minmod( W3, M3 );
	_W4 = minmod( W4, M4 );
	_W5 = minmod( W5, M5 );
	/* inverted transformation matrix [S-1] */
	/*_S21 =  w_hh;  _S31 =  u_hh; _S41 =  v_hh; */
	/*_S34 =  u_ee; _S44 =  w_ee; */
	_S24 = w_ee + ll; _S25 = w_ee - ll;
	_S54 = _S55 = nn + _1_2_K_1;
	_S54 += mm; _S55 -= mm;
	/* parameters' vector at the nearer interface */
	k1 = _M1 * C1 + _W1 * C2;
	k2 = _M2 * C1 + _W2 * C2;
	k3 = _M3 * C1 + _W3 * C2;
	k4 = _M4 * C1 + _W4 * C2;
	k5 = _M5 * C1 + _W5 * C2;
	kk = k4 + k5;
	zU1[_i][_j][ k] = u1 +   hh * k1 + 	ee  * kk;
	zU4[_i][_j][ k] = u4 + w_hh * k1 + 			   _S24 * k4 + _S25 * k5;
	zU2[_i][_j][ k] = u2 + u_hh * k1 + 	  k2 		 + u_ee * kk;
	zU3[_i][_j][ k] = u3 + v_hh * k1 + 				  k3 + w_ee * kk;
	zU5[_i][_j][ k] = u5 +   gg * k1 +    U * k2 +    V * k3 + _S54 * k4 + _S55 * k5;
	/* ... farther interface */
	k1 = _M1 * C2 + _W1 * C1;
	k2 = _M2 * C2 + _W2 * C1;
	k3 = _M3 * C2 + _W3 * C1;
	k4 = _M4 * C2 + _W4 * C1;
	k5 = _M5 * C2 + _W5 * C1;
	kk =  k4 + k5;
	U1z[_i][_j][_k] = u1 -   hh * k1 - ee * kk;
	U4z[_i][_j][_k] = u4 - w_hh * k1 -  _S24 * k4 - _S25 * k5;
	U2z[_i][_j][_k] = u2 - u_hh * k1 -  k2  - u_ee * kk;
	U3z[_i][_j][_k] = u3 - v_hh * k1 -  k3 - w_ee * kk;
	U5z[_i][_j][_k] = u5 -   gg * k1 -  U * k2 -  V * k3 - _S54 * k4 - _S55 * k5;
      } /* end for */
    } /* end for */
  } /* end for */
} /* end Reconstruction() */

/***********************
* BOUNCONDONINTERFACES *  Completing the "outer" parameters' values at the domain boundary
***********************/
void BounCondOnInterfaces(void)
{
  /* x */
  for( j = 0; j < HIG; j++ ) {
    for( k = 0; k < DEP; k++ ) {
      /* left side - PERIODIC INFLOW */
      xU1[0][j][k] = xU1[LEN][j][k];
      xU2[0][j][k] = xU2[LEN][j][k];
      xU3[0][j][k] = xU3[LEN][j][k];
      xU4[0][j][k] = xU4[LEN][j][k];
      xU5[0][j][k] = xU5[LEN][j][k];
      /* right side - PERIODIC OUTFLOW */
      U1x[LEN][j][k] = U1x[0][j][k];
      U2x[LEN][j][k] = U2x[0][j][k];
      U3x[LEN][j][k] = U3x[0][j][k];
      U4x[LEN][j][k] = U4x[0][j][k];
      U5x[LEN][j][k] = U5x[0][j][k];
    } /* end for */
  } /* end for */
  /* y */
  for( i = 0; i < LEN; i++ ) {
    for( k = 0; k < DEP; k++ ) {
      /* bottom side - NO SLIP */
      yU1[i][0][k] =   U1y[i][0][k];
      yU2[i][0][k] = - U2y[i][0][k];
      yU3[i][0][k] = - U3y[i][0][k];
      yU4[i][0][k] = - U4y[i][0][k];
      yU5[i][0][k] =   U5y[i][0][k];
      /* top side - NO SLIP */
      U1y[i][HIG][k] =   yU1[i][HIG][k];
      U2y[i][HIG][k] = - yU2[i][HIG][k];
      U3y[i][HIG][k] = - yU3[i][HIG][k];
      U4y[i][HIG][k] = - yU4[i][HIG][k];
      U5y[i][HIG][k] =   yU5[i][HIG][k];
    } /* end for */
  }  /* end for */
  /* z */
  for( i = 0; i < LEN; i++ ) {
    for( j = 0; j < HIG; j++ ) {
      /* back side - NO SLIP */
      zU1[i][j][0] =   U1z[i][j][0];
      zU2[i][j][0] = - U2z[i][j][0];
      zU3[i][j][0] = - U3z[i][j][0];
      zU4[i][j][0] = - U4z[i][j][0];
      zU5[i][j][0] =   U5z[i][j][0];
      /* front side - NO SLIP */
      U1z[i][j][DEP] =   zU1[i][j][DEP];
      U2z[i][j][DEP] = - zU2[i][j][DEP];
      U3z[i][j][DEP] = - zU3[i][j][DEP];
      U4z[i][j][DEP] = - zU4[i][j][DEP];
      U5z[i][j][DEP] =   zU5[i][j][DEP];
    } /* end for */
  } /* end for */

} /* end BounCondOnInterfaces( void ) */

/***********
*  FLUXES  *   Gasdynamical+molecular fluxes at the interfaces
***********/
void Fluxes( void )
{
/*register*/ unsigned i, j, k, i_, j_, k_, i__, j__, k__;
/*register*/ float RU; /* convective normal mass flux */

	/*--- X-fluxes ---*/
  for( i = 0, i_ = 1; i < LENN; i++, i_++ ) {
    for( j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++ ) {
      for( k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++ ) {
				/*--- Inviscid fluxes ---*/
	/* left parameters */
	_R =  xU1[i][j][k];
	_U =  xU2[i][j][k] / _R;
	_V =  xU3[i][j][k] / _R;
	_W =  xU4[i][j][k] / _R;
	_P = (xU5[i][j][k] - 0.5 * _R * (_U*_U + _V*_V + _W*_W)) * K_1;
	_C  = sqrt(K * _P / _R);
	_jo = _K_1 * _P; /* _P - _C * _C * _R */
	_jp = _U + ( RU = _P / ( _R * _C ) );
	_jm = _U -   RU;
	/* right parameters */
	R_ =   U1x[i][j][k];
	U_ =   U2x[i][j][k] / R_;
	V_ =   U3x[i][j][k] / R_;
	W_ =   U4x[i][j][k] / R_;
	P_ = ( U5x[i][j][k]
	       - 0.5 * R_ * ( U_*U_ + V_*V_ + W_*W_ )
	       ) * K_1;
	C_  = sqrt( K * P_ / R_ );
	jo_ = _K_1 * P_; /* P_ - C_ * C_ * R_ */
	jp_ = U_ + ( RU = P_ / ( R_ * C_ ) );
	jm_ = U_ -   RU;
	/* linearized procedure */
	U = 0.5 * ( _U + U_ );
	C = 0.5 * ( _C + C_ );
	C_p = U + C; C_m = U - C; C_o = U;
	/*              labelx:
	 */
	/* J_p and al_p */
	if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
	else           { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
	/* J_m and al_m */
	if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
	else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
	/* J_o and al_o */
	if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
	else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
	/* parameters at the interface */
	P = ( J_p - J_m ) / ( al_p - al_m );
	U = J_p - al_p * P;
	R = ( J_o - P ) / al_o;
	/* check signs of the Ch-velocities (optional) */
	/* C = sqrt( K * P / R );
	   C_p_ = U + C; C_m_ = U - C; C_o_ = U;
	   if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
	   C_p = C_p_; C_m = C_m_; C_o = C_o_;
	   goto labelx;
	   }
	*/
	/* tangential velocities */
	if( U > 0 ) { V = _V; W = _W; }
	else {        V = V_; W = W_; }
				/*--- Viscous fluxes ---*/
	/*-- forward cells --*/
	/* forward */
	R_ =    U1_[i_ ][j_ ][k_ ];
	U_ =    U2_[i_ ][j_ ][k_ ] / R_;
	V_ =    U3_[i_ ][j_ ][k_ ] / R_;
	W_ =    U4_[i_ ][j_ ][k_ ] / R_;
	P_ = (  U5_[i_ ][j_ ][k_ ]
		- 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ )
		) * K_1;
	T_ = P_ / ( R_VOZD * R_ );
	/* forward upper */
	RU = 1./U1_[i_ ][j__][k_ ];
	Uj =    U2_[i_ ][j__][k_ ] * RU;
	Vj =    U3_[i_ ][j__][k_ ] * RU;
	Wj =    U4_[i_ ][j__][k_ ] * RU;
	/* forward lower */
	RU = 1./U1_[i_ ][j  ][k_ ];
	jU =    U2_[i_ ][j  ][k_ ] * RU;
	jV =    U3_[i_ ][j  ][k_ ] * RU;
	jW =    U4_[i_ ][j  ][k_ ] * RU;
	/* forward nearer */
	RU = 1./U1_[i_ ][j_ ][k__];
	Uk =    U2_[i_ ][j_ ][k__] * RU;
	Vk =    U3_[i_ ][j_ ][k__] * RU;
	Wk =    U4_[i_ ][j_ ][k__] * RU;
	/* forward farther */
	RU = 1./U1_[i_ ][j_ ][k  ];
	kU =    U2_[i_ ][j_ ][k  ] * RU;
	kV =    U3_[i_ ][j_ ][k  ] * RU;
	kW =    U4_[i_ ][j_ ][k  ] * RU;
	/*-- backward cells -- */
	/* backward */
	_R =    U1_[i  ][j_ ][k_ ];
	_U =    U2_[i  ][j_ ][k_ ] / _R;
	_V =    U3_[i  ][j_ ][k_ ] / _R;
	_W =    U4_[i  ][j_ ][k_ ] / _R;
	_P = (  U5_[i  ][j_ ][k_ ]
		- 0.5 * _R * ( _U * _U + _V * _V + _W * _W )
		) * K_1;
	_T = _P / ( R_VOZD * _R );
	/* backward upper */
	RU = 1./U1_[i  ][j__][k_ ];
	uj =    U2_[i  ][j__][k_ ] * RU;
	vj =    U3_[i  ][j__][k_ ] * RU;
	wj =    U4_[i  ][j__][k_ ] * RU;
	/* backward lower */
	RU = 1./U1_[i  ][j  ][k_ ];
	ju =    U2_[i  ][j  ][k_ ] * RU;
	jv =    U3_[i  ][j  ][k_ ] * RU;
	jw =    U4_[i  ][j  ][k_ ] * RU;
	/* backward nearer */
	RU = 1./U1_[i  ][j_ ][k__];
	uk =    U2_[i  ][j_ ][k__] * RU;
	vk =    U3_[i  ][j_ ][k__] * RU;
	wk =    U4_[i  ][j_ ][k__] * RU;
	/* backward farther */
	RU = 1./U1_[i  ][j_ ][k  ];
	ku =    U2_[i  ][j_ ][k  ] * RU;
	kv =    U3_[i  ][j_ ][k  ] * RU;
	kw =    U4_[i  ][j_ ][k  ] * RU;
	/* derivatives of velocities */
	/* x */
	du_dx = ( U_ - _U ) * _deltaX;
	dv_dx = ( V_ - _V ) * _deltaX;
	dw_dx = ( W_ - _W ) * _deltaX;
	/* y */
	du_dy = ( Uj + uj - jU - ju ) * _4deltaY;
	dv_dy = ( Vj + vj - jV - jv ) * _4deltaY;
	dw_dy = ( Wj + wj - jW - jw ) * _4deltaY;
	/* z */
	du_dz = ( Uk + uk - kU - ku ) * _4deltaZ;
	dv_dz = ( Vk + vk - kV - kv ) * _4deltaZ;
	dw_dz = ( Wk + wk - kW - kw ) * _4deltaZ;
	/* mean velocities and density */
	rr = 0.5 * ( R_ + _R );
	uu = 0.5 * ( U_ + _U );
	vv = 0.5 * ( V_ + _V );
	ww = 0.5 * ( W_ + _W );
	/* stresses */
	S12 = du_dy + dv_dx;
	S13 = du_dz + dw_dx;
	S23 = dv_dz + dw_dy;
	_S_ = sqrt( 2. * ( du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
		    +  ( S12   * S12   + S13   * S13   + S23   * S23 ) );
	mu_E = mu_L + ( mu_T = rr * CsDDx[j][k] * _S_ ); /* mu_T = r * Cs * delta * delta * | S | */
	sigma_xx = 0.66666667 * mu_E * ( du_dx + du_dx - dv_dy - dw_dz );
	sigma_xy = mu_E * ( du_dy + dv_dx );
	sigma_xz = mu_E * ( du_dz + dw_dx );
	/* heat flux */
	lambda_E = lambda_L + mu_T * cp_Pr_T;
	q_x = - lambda_E * ( T_ - _T ) * _deltaX;
	/* summary X-fluxes evaluation */
	Fx1[i][j][k] = RU = R * U;
	Fx2[i][j][k] = RU * U + P - sigma_xx;
	Fx3[i][j][k] = RU * V     - sigma_xy;
	Fx4[i][j][k] = RU * W	  - sigma_xz;
	Fx5[i][j][k] = U * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) )
	  - uu * sigma_xx - vv * sigma_xy - ww * sigma_xz + q_x;
      }
    }
  }

  /*--- Y-fluxes ---*/
  for( j = 0, j_ = 1; j < HIGG; j++, j_++ ) {
    for( i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++ ) {
      for( k = 0, k_ = 1, k__ = 2; k < DEP; k++, k_++, k__++ ) {
	/*--- Inviscid fluxes ---*/
	/* lower */
	_R =   yU1[i][j][k];
	_U =   yU2[i][j][k] / _R;
	_V =   yU3[i][j][k] / _R;
	_W =   yU4[i][j][k] / _R;
	_P = ( yU5[i][j][k]
	       - 0.5 * _R * ( _U*_U + _V*_V + _W*_W )
	       ) * K_1;
	_C  = sqrt( K * _P / _R );
	_jo = _K_1 * _P; /* _P - _C * _C * _R */
	_jp = _V + ( RU = _P / ( _R * _C ) );
	_jm = _V -   RU;
	/* upper */
	R_ =   U1y[i][j][k];
	U_ =   U2y[i][j][k] / R_;
	V_ =   U3y[i][j][k] / R_;
	W_ =   U4y[i][j][k] / R_;
	P_ = ( U5y[i][j][k]
	       - 0.5 * R_ * ( U_*U_ + V_*V_ + W_*W_ )
	       ) * K_1;
	C_  = sqrt( K * P_ / R_ );
	jo_ = _K_1 * P_; /* P_ - C_ * C_ * R_ */
	jp_ = V_ + ( RU = P_ / ( R_ * C_ ) );
	jm_ = V_ -   RU;
	/* linearized procedure */
	V = 0.5 * ( _V + V_ );
	C = 0.5 * ( _C + C_ );
	C_p = V + C; C_m = V - C; C_o = V;
	/*              labely:
	 */
	/* Jp and al_p */
	if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
	else           { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
	/* Jm and al_m */
	if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
	else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
	/* Jo and al_o */
	if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
	else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
	/* parameters at the interface */
	P = ( J_p - J_m ) / ( al_p - al_m );
	V = J_p - al_p * P;
	R = ( J_o - P ) / al_o;
	/* check signs of the Ch-velocities */
	/*
	  C = sqrt( K * P / R );
	  C_p_ = V + C; C_m_ = V - C; C_o_ = V;
	  if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
	  C_p = C_p_; C_m = C_m_; C_o = C_o_;
	  goto labely;
	  }
	*/
	/* tangential velocities */
	if( V > 0 ) { U = _U; W = _W; }
	else {        U = U_; W = W_; }
				/*--- Viscous fluxes ---*/
	/*-- upper cells --*/
	/* upper */
	R_ =    U1_[i_ ][j_ ][k_ ];
	U_ =    U2_[i_ ][j_ ][k_ ] / R_;
	V_ =    U3_[i_ ][j_ ][k_ ] / R_;
	W_ =    U4_[i_ ][j_ ][k_ ] / R_;
	P_ = (  U5_[i_ ][j_ ][k_ ]
		- 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ )
		) * K_1;
	T_ = P_ / ( R_VOZD * R_ );
	/* upper nearer */
	RU = 1./U1_[i_ ][j_ ][k__];
	Uk =    U2_[i_ ][j_ ][k__] * RU;
	Vk =    U3_[i_ ][j_ ][k__] * RU;
	Wk =    U4_[i_ ][j_ ][k__] * RU;
	/* upper farther */
	RU = 1./U1_[i_ ][j_ ][k  ];
	kU =    U2_[i_ ][j_ ][k  ] * RU;
	kV =    U3_[i_ ][j_ ][k  ] * RU;
	kW =    U4_[i_ ][j_ ][k  ] * RU;
	/* upper forward */
	RU = 1./U1_[i__][j_ ][k_ ];
	Ui =    U2_[i__][j_ ][k_ ] * RU;
	Vi =    U3_[i__][j_ ][k_ ] * RU;
	Wi =    U4_[i__][j_ ][k_ ] * RU;
	/* upper backward */
	RU = 1./U1_[i  ][j_ ][k_ ];
	iU =    U2_[i  ][j_ ][k_ ] * RU;
	iV =    U3_[i  ][j_ ][k_ ] * RU;
	iW =    U4_[i  ][j_ ][k_ ] * RU;
	/*-- lower cells --*/
	/* lower */
	_R =    U1_[i_ ][j  ][k_ ];
	_U =    U2_[i_ ][j  ][k_ ] / _R;
	_V =    U3_[i_ ][j  ][k_ ] / _R;
	_W =    U4_[i_ ][j  ][k_ ] / _R;
	_P = (  U5_[i_ ][j  ][k_ ]
		- 0.5 * _R * ( _U * _U + _V * _V + _W * _W )
		) * K_1;
	_T = _P / ( R_VOZD * _R );
	/* lower nearer */
	RU = 1./U1_[i_ ][j  ][k__];
	uk =    U2_[i_ ][j  ][k__] * RU;
	vk =    U3_[i_ ][j  ][k__] * RU;
	wk =    U4_[i_ ][j  ][k__] * RU;
	/* lower farther */
	RU = 1./U1_[i_ ][j  ][k  ];
	ku = 	U2_[i_ ][j  ][k  ] * RU;
	kv = 	U3_[i_ ][j  ][k  ] * RU;
	kw = 	U4_[i_ ][j  ][k  ] * RU;
	/* lower forward */
	RU = 1./U1_[i__][j  ][k_ ];
	ui =    U2_[i__][j  ][k_ ] * RU;
	vi =    U3_[i__][j  ][k_ ] * RU;
	wi =    U4_[i__][j  ][k_ ] * RU;
	/* lower backward */
	RU = 1./U1_[i  ][j  ][k_ ];
	iu =    U2_[i  ][j  ][k_ ] * RU;
	iv =    U3_[i  ][j  ][k_ ] * RU;
	iw =    U4_[i  ][j  ][k_ ] * RU;
	/* derivatives of velocities */
	/* x */
	du_dx = ( Ui + ui - iU - iu ) * _4deltaX;
	dv_dx = ( Vi + vi - iV - iv ) * _4deltaX;
	dw_dx = ( Wi + wi - iW - iw ) * _4deltaX;
	/* y */
	du_dy = ( U_ - _U ) * _deltaY;
	dv_dy = ( V_ - _V ) * _deltaY;
	dw_dy = ( W_ - _W ) * _deltaY;
	/* z */
	du_dz = ( Uk + uk - kU - ku ) * _4deltaZ;
	dv_dz = ( Vk + vk - kV - kv ) * _4deltaZ;
	dw_dz = ( Wk + wk - kW - kw ) * _4deltaZ;
	/* mean velocities and density */
	rr = 0.5 * ( R_ + _R );
	uu = 0.5 * ( U_ + _U );
	vv = 0.5 * ( V_ + _V );
	ww = 0.5 * ( W_ + _W );
	/* stresses */
	S12 = du_dy + dv_dx;
	S13 = du_dz + dw_dx;
	S23 = dv_dz + dw_dy;
	_S_ = sqrt( 2. * ( du_dx*du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
		    +  ( S12  * S12  + S13   * S13   + S23   * S23 ) );
	mu_E = mu_L + ( mu_T = rr * CsDDy[j][k] * _S_ );  /* mu_T = r * Cs * delta * delta * | S |  */
	sigma_yx = mu_E * ( dv_dx + du_dy );
	sigma_yy = 0.66666667 * mu_E * ( dv_dy + dv_dy - du_dx - dw_dz );
	sigma_yz = mu_E * ( dv_dz + dw_dy );
	/* heat flux */
	lambda_E = lambda_L + mu_T * cp_Pr_T;
	q_y = - lambda_E * ( T_ - _T ) * _deltaY;
	/* summary Y-fluxes evaluation */
	Fy1[i][j][k] = RU = R * V;
	Fy2[i][j][k] = RU * U     - sigma_yx;
	Fy3[i][j][k] = RU * V + P - sigma_yy;
	Fy4[i][j][k] = RU * W     - sigma_yz;
	Fy5[i][j][k] = V * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) )
	  - uu * sigma_yx - vv * sigma_yy - ww * sigma_yz + q_y;
      }
    }
  }

  /*--- Z-fluxes ---*/
  for (k = 0, k_ = 1; k < DEPP; k++, k_++) {
    for (j = 0, j_ = 1, j__ = 2; j < HIG; j++, j_++, j__++) {
      for (i = 0, i_ = 1, i__ = 2; i < LEN; i++, i_++, i__++) {
	/*--- Inviscid fluxes ---*/
	/* farther */
	_R =   zU1[i][j][k];
	_U =   zU2[i][j][k] / _R;
	_V =   zU3[i][j][k] / _R;
	_W =   zU4[i][j][k] / _R;
	_P = ( zU5[i][j][k]
	       - 0.5 * _R * ( _U*_U + _V*_V + _W*_W )
	       ) * K_1;
	_C  = sqrt( K * _P / _R );
	_jo = _K_1 * _P;                      /* _P - _C * _C * _R; */
	_jp = _W + ( RU = _P / ( _R * _C ) );
	_jm = _W -   RU;
	/* nearer */
	R_ =  U1z[i][j][k];
	U_ =  U2z[i][j][k] / R_;
	V_ =  U3z[i][j][k] / R_;
	W_ =  U4z[i][j][k] / R_;
	P_ = (U5z[i][j][k] - 0.5 * R_ * (U_*U_ + V_*V_ + W_*W_)) * K_1;
	C_  = sqrt( K * P_ / R_ );
	jo_ = _K_1 * P_; 			           /* P_ - C_ * C_ * R_; */
	jp_ = W_ + ( RU = P_ / ( R_ * C_ ) );
	jm_ = W_ -   RU;
	/* linearized procedure */
	W = 0.5 * (_W + W_);
	C = 0.5 * (_C + C_);
	C_p = W + C; C_m = W - C; C_o = W;
	/*              labelz:
	 */
	/* Jp and al_p */
	if( C_p > 0. ) { al_p = K_1_2 * ( _jm - _jp ) / _jo; J_p  = _jp; }
	else 	       { al_p = K_1_2 * ( jm_ - jp_ ) / jo_; J_p  = jp_; }
	/* Jm and al_m */
	if( C_m > 0. ) { al_m = K_1_2 * ( _jp - _jm ) / _jo; J_m  = _jm; }
	else           { al_m = K_1_2 * ( jp_ - jm_ ) / jo_; J_m  = jm_; }
	/* Jo and al_o */
	if( C_o > 0. ) { al_o = _05K * ( _jp - _jm ); al_o = - al_o * al_o; J_o  = _jo; }
	else           { al_o = _05K * ( jp_ - jm_ ); al_o = - al_o * al_o; J_o  = jo_; }
	/* parameters at the interface */
	P = ( J_p - J_m ) / ( al_p - al_m );
	W = J_p - al_p * P;
	R = ( J_o - P ) / al_o;
	/* check signs of the Ch-velocities */
	/*
	  C = sqrt( K * P / R );
	  C_p_ = W + C; C_m_ = W - C; C_o_ = W;
	  if( C_p_*C_p < 0 || C_m_*C_m < 0 || C_o_*C_o < 0 ) {
	  C_p = C_p_; C_m = C_m_; C_o = C_o_;
	  goto labelz;
	  }
	*/
	/* tangential velocities */
	if( W > 0 ) { U = _U; V = _V; }
	else {        U = U_; V = V_; }
				/*--- Viscous fluxes ---*/
	/*-- nearer cells --*/
	/* nearer */
	R_ =    U1_[i_ ][j_ ][k_ ];
	U_ =    U2_[i_ ][j_ ][k_ ] / R_;
	V_ =    U3_[i_ ][j_ ][k_ ] / R_;
	W_ =    U4_[i_ ][j_ ][k_ ] / R_;
	P_ = (  U5_[i_ ][j_ ][k_ ]
		- 0.5 * R_ * ( U_ * U_ + V_ * V_ + W_ * W_ )
		) * K_1;
	T_ = P_ / ( R_VOZD * R_ );
	/* nearer forward */
	RU = 1./U1_[i__][j_ ][k_ ];
	Ui =    U2_[i__][j_ ][k_ ] * RU;
	Vi =    U3_[i__][j_ ][k_ ] * RU;
	Wi =    U4_[i__][j_ ][k_ ] * RU;
	/* nearer backward */
	RU = 1./U1_[i  ][j_ ][k_ ];
	iU =    U2_[i  ][j_ ][k_ ] * RU;
	iV =    U3_[i  ][j_ ][k_ ] * RU;
	iW =    U4_[i  ][j_ ][k_ ] * RU;
	/* nearer upper */
	RU = 1./U1_[i_ ][j__][k_ ];
	Uj =    U2_[i_ ][j__][k_ ] * RU;
	Vj =    U3_[i_ ][j__][k_ ] * RU;
	Wj =    U4_[i_ ][j__][k_ ] * RU;
	/* nearer lower */
	RU = 1./U1_[i_ ][j  ][k_ ];
	jU =    U2_[i_ ][j  ][k_ ] * RU;
	jV =    U3_[i_ ][j  ][k_ ] * RU;
	jW =    U4_[i_ ][j  ][k_ ] * RU;
	/*-- farther cells --*/
	/* farther */
	_R =    U1_[i_ ][j_ ][k  ];
	_U =    U2_[i_ ][j_ ][k  ] / _R;
	_V =    U3_[i_ ][j_ ][k  ] / _R;
	_W =    U4_[i_ ][j_ ][k  ] / _R;
	_P = (  U5_[i_ ][j_ ][k  ]
		- 0.5 * _R * ( _U * _U + _V * _V + _W * _W )
		) * K_1;
	_T = _P / ( R_VOZD * _R );
	/* farther forward */
	RU = 1./U1_[i__][j_ ][k  ];
	ui =    U2_[i__][j_ ][k  ] * RU;
	vi =    U3_[i__][j_ ][k  ] * RU;
	wi =    U4_[i__][j_ ][k  ] * RU;
	/* farther backward */
	RU = 1./U1_[i  ][j_ ][k  ];
	iu =    U2_[i  ][j_ ][k  ] * RU;
	iv =    U3_[i  ][j_ ][k  ] * RU;
	iw =    U4_[i  ][j_ ][k  ] * RU;
	/* farther upper */
	RU = 1./U1_[i_ ][j__][k  ];
	uj =    U2_[i_ ][j__][k  ] * RU;
	vj =    U3_[i_ ][j__][k  ] * RU;
	wj =    U4_[i_ ][j__][k  ] * RU;
	/* farther lower */
	RU = 1./U1_[i_ ][j  ][k  ];
	ju =    U2_[i_ ][j  ][k  ] * RU;
	jv =    U3_[i_ ][j  ][k  ] * RU;
	jw =    U4_[i_ ][j  ][k  ] * RU;
	/* derivatives of velocities */
	/* x */
	du_dx = ( Ui + ui - iU - iu ) * _4deltaX;
	dv_dx = ( Vi + vi - iV - iv ) * _4deltaX;
	dw_dx = ( Wi + wi - iW - iw ) * _4deltaX;
	/* y */
	du_dy = ( Uj + uj - jU - ju ) * _4deltaY;
	dv_dy = ( Vj + vj - jV - jv ) * _4deltaY;
	dw_dy = ( Wj + wj - jW - jw ) * _4deltaY;
	/* z */
	du_dz = ( U_ - _U ) * _deltaZ;
	dv_dz = ( V_ - _V ) * _deltaZ;
	dw_dz = ( W_ - _W ) * _deltaZ;
	/* mean velocities and density */
	rr = 0.5 * ( R_ + _R );
	uu = 0.5 * ( U_ + _U );
	vv = 0.5 * ( V_ + _V );
	ww = 0.5 * ( W_ + _W );
	/* stresses */
	S12 = du_dy + dv_dx;
	S13 = du_dz + dw_dx;
	S23 = dv_dz + dw_dy;
	_S_ = sqrt( 2. * ( du_dx * du_dx + dv_dy * dv_dy + dw_dz * dw_dz )
		    +  ( S12   * S12   + S13   * S13   + S23   * S23 ) );
	mu_E = mu_L + ( mu_T = rr * CsDDz[j][k] * _S_ ); /* mu_T = r * Cs * delta * delta * | S | */
	sigma_zx = mu_E * ( dw_dx + du_dz );
	sigma_zy = mu_E * ( dw_dy + dv_dz );
	sigma_zz = 0.66666667 * mu_E * ( dw_dz + dw_dz - du_dx - dv_dy );
	/* heat flux */
	lambda_E = lambda_L + mu_T * cp_Pr_T;
	q_z = - lambda_E * ( T_ - _T ) * _deltaZ;
	/* summary Z-fluxes evaluation */
	Fz1[i][j][k] = RU = R * W;
	Fz2[i][j][k] = RU * U 	  - sigma_zx;
	Fz3[i][j][k] = RU * V 	  - sigma_zy;
	Fz4[i][j][k] = RU * W + P - sigma_zz;
	Fz5[i][j][k] = W * ( K_K_1 * P + 0.5 * R * ( U * U + V * V + W * W ) )
	  - uu * sigma_zx - vv * sigma_zy - ww * sigma_zz + q_z;
      }
    }
  }

} /* end Fluxes() */

/**************
*  EVOLUTION  *    Compute new/(intermediate) values of conservative variables
**************/
void Evolution(void)
{
/*register*/ unsigned i, j, k, _i, _j, _k; /* i_, j_, k_*/
/*register*/ float R, RU;

  if (Stage == 1) {
    for (i = 1, _i = 0; i < LENN; i++, _i++) {
      for (j = 1, _j = 0; j < HIGG; j++, _j++) { /* indices i-1, k-1 etc */
	for (k = 1, _k = 0; k < DEPP; k++, _k++) {
	  U1p[i][j][k] = ( R = U1[i][j][k] )
	    + deltaT_X * ( Fx1[_i][_j][_k] - Fx1[i][_j][_k] )
	    + deltaT_Y * ( Fy1[_i][_j][_k] - Fy1[_i][j][_k] )
	    + deltaT_Z * ( Fz1[_i][_j][_k] - Fz1[_i][_j][k] );
	  U2p[i][j][k] = ( RU = U2[i][j][k] )
	    + deltaT_X * ( Fx2[_i][_j][_k] - Fx2[i][_j][_k] )
	    + deltaT_Y * ( Fy2[_i][_j][_k] - Fy2[_i][j][_k] )
	    + deltaT_Z * ( Fz2[_i][_j][_k] - Fz2[_i][_j][k] )
	    + R * gdeltaT;
	  U3p[i][j][k] = U3[i][j][k]
	    + deltaT_X * ( Fx3[_i][_j][_k] - Fx3[i][_j][_k] )
	    + deltaT_Y * ( Fy3[_i][_j][_k] - Fy3[_i][j][_k] )
	    + deltaT_Z * ( Fz3[_i][_j][_k] - Fz3[_i][_j][k] );
	  U4p[i][j][k] = U4[i][j][k]
	    + deltaT_X * ( Fx4[_i][_j][_k] - Fx4[i][_j][_k] )
	    + deltaT_Y * ( Fy4[_i][_j][_k] - Fy4[_i][j][_k] )
	    + deltaT_Z * ( Fz4[_i][_j][_k] - Fz4[_i][_j][k] );
	  U5p[i][j][k] = U5[i][j][k]
	    + deltaT_X * ( Fx5[_i][_j][_k] - Fx5[i][_j][_k] )
	    + deltaT_Y * ( Fy5[_i][_j][_k] - Fy5[_i][j][_k] )
	    + deltaT_Z * ( Fz5[_i][_j][_k] - Fz5[_i][_j][k] )
	    + RU * gdeltaT;
	}
      }
    }
    /* next stage */
    Stage = 2;
    U1_ = U1p, U2_ = U2p, U3_ = U3p, U4_ = U4p, U5_ = U5p;
  } /* end if() First stage */

  else {
    for (i = 1, _i = 0; i < LENN; i++, _i++) {
      for (j = 1, _j = 0; j < HIGG; j++, _j++) { /* indices i-1, k-1 etc */
	for (k = 1, _k = 0; k < DEPP; k++, _k++) {
	  U1[i][j][k] = 0.5 * ( U1[i][j][k] + ( R = U1p[i][j][k] )
				+ deltaT_X * ( Fx1[_i][_j][_k] - Fx1[i][_j][_k] )
				+ deltaT_Y * ( Fy1[_i][_j][_k] - Fy1[_i][j][_k] )
				+ deltaT_Z * ( Fz1[_i][_j][_k] - Fz1[_i][_j][k] ) );
	  U2[i][j][k] = 0.5 * ( U2[i][j][k] + ( RU = U2p[i][j][k] )
				+ deltaT_X * ( Fx2[_i][_j][_k] - Fx2[i][_j][_k] )
				+ deltaT_Y * ( Fy2[_i][_j][_k] - Fy2[_i][j][_k] )
				+ deltaT_Z * ( Fz2[_i][_j][_k] - Fz2[_i][_j][k] )
				+ R * gdeltaT  );
	  U3[i][j][k] = 0.5 * ( U3[i][j][k] + U3p[i][j][k]
				+ deltaT_X * ( Fx3[_i][_j][_k] - Fx3[i][_j][_k] )
				+ deltaT_Y * ( Fy3[_i][_j][_k] - Fy3[_i][j][_k] )
				+ deltaT_Z * ( Fz3[_i][_j][_k] - Fz3[_i][_j][k] ) );
	  U4[i][j][k] = 0.5 * ( U4[i][j][k] + U4p[i][j][k]
				+ deltaT_X * ( Fx4[_i][_j][_k] - Fx4[i][_j][_k] )
				+ deltaT_Y * ( Fy4[_i][_j][_k] - Fy4[_i][j][_k] )
				+ deltaT_Z * ( Fz4[_i][_j][_k] - Fz4[_i][_j][k] ) );
	  U5[i][j][k] = 0.5 * ( U5[i][j][k] + U5p[i][j][k]
				+ deltaT_X * ( Fx5[_i][_j][_k] - Fx5[i][_j][_k] )
				+ deltaT_Y * ( Fy5[_i][_j][_k] - Fy5[_i][j][_k] )
				+ deltaT_Z * ( Fz5[_i][_j][_k] - Fz5[_i][_j][k] )
				+ RU * gdeltaT );
	}
      }
    }

    Stage = 1;
    U1_ = U1, U2_ = U2, U3_ = U3, U4_ = U4, U5_ = U5;
  } /* end else Second stage */

} /* end Evolution() */


/*
  Write readings of velocity probes to the binary "probes.dat" file
  and cross-sectional fluxes and wall friction stress to ASCII "probes.txt"
*/

/***********
*  PROBES  *
***********/
void Probes(void)
{

  if (Stage == 1) {

    /*--- mean fluxes at this stage ---*/
    i = LEN/2;
    for (j = 0, M = I = E = 0.; j < HIG; j++) {
      for (k = 0; k < DEP; k++) {
	M += Fx1[i][j][k]; /* mass     */
	I += Fx2[i][j][k]; /* momentum */
	E += Fx5[i][j][k]; /* energy   */
      } /* end for */
    } /* end for */

    /*--- mean wall stresses at this stage ---*/
    tauZ = tauY = 0.;
    /* wall sigma_zx stresses */
    for (j = 0; j < HIG; j++)
      tauZ += Fz2[i][j][0] - Fz2[i][j][DEP];
    /* wall sigma_yx stresses */
    for (k = 0; k < DEP; k++)
      tauY += Fy2[i][0][k] - Fy2[i][HIG][k];
  }

  else {

    /*--- plus mean fuxes at this stage ---*/
    i = LEN/2;
    for (j = 0; j < HIG; j++) {
      for (k = 0; k < DEP; k++) {
	M += Fx1[i][j][k]; /* mass     */
	I += Fx2[i][j][k]; /* momentum */
	E += Fx5[i][j][k]; /* energy   */
      } /* end for */
    } /* end for */

    /*--- plus mean wall stresses at this stage ---*/
    /* wall sigma_zx stresses */
    for (j = 0; j < HIG; j++)
      tauZ += Fz2[i][j][0] - Fz2[i][j][DEP];
    /* wall sigma_yx stresses */
    for (k = 0; k < DEP; k++)
      tauY += Fy2[i][0][k] - Fy2[i][HIG][k];

    /*--- write mean fluxes and mean wall stess ---*/
    fprintf(pFprobest, "%d %.5f %.1f %.0f %f\n", step,
	    0.5*M/(DEP*HIG), 0.5*I/(DEP*HIG), 0.5*E/(DEP*HIG),
	    0.5*(tauZ*deltaY + tauZ*deltaZ)
	    / (2. * (DEP*deltaZ + HIG*deltaY)));

    /*--- number of current step ---*/
    fwrite(&step, sizeof(unsigned long), 1, pFprobesb);

    /*--- writes HIG-sized blocks of current R and RU ---*/
    for (j = 1, k = DEP/2, i = LEN/2; j < HIGG; j++)
      probes[j-1] = U1[i][j][k];
    fwrite(probes, sizeof(float)*HIG, 1, pFprobesb);
    for (j = 1; j < HIGG; j++)
      probes[j-1] = U2[i][j][k];
    fwrite(probes, sizeof(float)*HIG, 1, pFprobesb);

  } /* end else */

} /* end Probes() */


/***********
*  OUTPUT  *   Outputs vorticity field to ASCII file
***********/
void Output(void)
{
  sprintf(str, "%d", step);
  if ((pF = fopen(str, "w")) == NULL) {
    printf("Cannot open file \"%s\" ", str);
    exit(-1);
  }

  //  fprintf(pF, "%d\n%d\n%f\n%f\n", LEN, HIG, deltaX, deltaY);
  k = DEP/2;
  for (i = 1; i < LENN; i++) {
    for (j = 1; j < HIGG; j++) {
      dv_dx = (U3[i+1][j][k]/U1[i+1][j][k] - U3[i-1][j][k]/U1[i-1][j][k]) * _2deltaX;
      dw_dx = (U4[i+1][j][k]/U1[i+1][j][k] - U4[i-1][j][k]/U1[i-1][j][k]) * _2deltaX;
      du_dy = (U2[i][j+1][k]/U1[i][j+1][k] - U2[i][j-1][k]/U1[i][j-1][k]) * _2deltaY;
      dw_dy = (U4[i][j+1][k]/U1[i][j+1][k] - U4[i][j-1][k]/U1[i][j-1][k]) * _2deltaY;
      du_dz = (U2[i][j][k+1]/U1[i][j][k+1] - U2[i][j][k-1]/U1[i][j][k-1]) * _2deltaZ;
      dv_dz = (U3[i][j][k+1]/U1[i][j][k+1] - U3[i][j][k-1]/U1[i][j][k-1]) * _2deltaZ;
      U = 0.5 * (dw_dy - dv_dz);
      V = 0.5 * (du_dz - dw_dx);
      W = 0.5 * (dv_dx - du_dy);
      fprintf( pF, "%f \n", sqrt( U * U + V * V + W * W ) );
    } /* end for */
    fprintf( pF, "\n");
  } /* end for */
  fclose( pF );

} /* end Output() */

/***************
*  OUTPUT_VTK  *   Outputs vorticity field to ASCII VTK file
***************/
/** Truong added 26/02/2018 **/
void Output_VTK(void)
{

  float coordx[LENN];
  float coordy[HIGG];

  for (int ii = 0; ii < LENN; ii++) {
    coordx[ii] = ii * deltaX;
  }


  for (int jj = 0; jj < HIGG; jj++) {
    coordy[jj] = (jj) * deltaY;
  }

  sprintf(str, "%d", step);
  char str2[20] = ".vtk";
  sprintf(str2,strcat(str,str2));
  if ((pF = fopen(str2, "w")) == NULL) {
    printf("Cannot open file \"%s\" ", str2);
    exit(-1);
  }

  //  fprintf(pF, "%d\n%d\n%f\n%f\n", LEN, HIG, deltaX, deltaY);
  fprintf(pF,"# vtk DataFile Version 2.0 \n");
  fprintf(pF,"LES slice \n");
  fprintf(pF,"ASCII \n");
  fprintf(pF,"DATASET STRUCTURED_GRID \n");

  fprintf(pF,"DIMENSIONS %d %d %d \n",LENN,HIGG,1 );
  fprintf(pF,"POINTS %d float \n",LENN*HIGG);


  for (int jj = 0; jj < HIGG; jj++)
  {
      for (int ii = 0; ii < LENN; ii++)
      {
          fprintf(pF,"%f %f %f \n",coordx[ii],coordy[jj],0.0);
      }
  }

  fprintf(pF,"CELL_DATA %d \n",(LENN-1)*(HIGG-1));
  fprintf(pF,"SCALARS vorticity float \n");
  fprintf(pF,"LOOKUP_TABLE default \n");

  k = DEP/2;

  for (i = 1; i < LENN; i++) {
    for (j = 1; j < HIGG; j++) {
      dv_dx = (U3[i+1][j][k]/U1[i+1][j][k] - U3[i-1][j][k]/U1[i-1][j][k]) * _2deltaX;
      dw_dx = (U4[i+1][j][k]/U1[i+1][j][k] - U4[i-1][j][k]/U1[i-1][j][k]) * _2deltaX;
      du_dy = (U2[i][j+1][k]/U1[i][j+1][k] - U2[i][j-1][k]/U1[i][j-1][k]) * _2deltaY;
      dw_dy = (U4[i][j+1][k]/U1[i][j+1][k] - U4[i][j-1][k]/U1[i][j-1][k]) * _2deltaY;
      du_dz = (U2[i][j][k+1]/U1[i][j][k+1] - U2[i][j][k-1]/U1[i][j][k-1]) * _2deltaZ;
      dv_dz = (U3[i][j][k+1]/U1[i][j][k+1] - U3[i][j][k-1]/U1[i][j][k-1]) * _2deltaZ;
      U = 0.5 * (dw_dy - dv_dz);
      V = 0.5 * (du_dz - dw_dx);
      W = 0.5 * (dv_dx - du_dy);
      fprintf( pF, "%f \n", sqrt( U * U + V * V + W * W ) );
    } /* end for */

  } /* end for */
  fclose( pF );

} /* end Output_VTK() */



/*************
*  Finalize  *     Saves solution arrays to the binary "backup.dat" file
*************/
void Finalize( void )
{
  if ((pF = fopen("backup.dat", "wb")) == NULL) {
    fprintf(stderr, "can't open \"backup.dat\" to write.\n");
    exit(-1);
  }
  for (i = 1; i < LENN; i++) {
    for (j = 1; j < HIGG; j++) {
      fwrite(&U1[i][j][1], DEP * sizeof(U1[i][j][1]), 1, pF);
      fwrite(&U2[i][j][1], DEP * sizeof(U2[i][j][1]), 1, pF);
      fwrite(&U3[i][j][1], DEP * sizeof(U3[i][j][1]), 1, pF);
      fwrite(&U4[i][j][1], DEP * sizeof(U4[i][j][1]), 1, pF);
      fwrite(&U5[i][j][1], DEP * sizeof(U5[i][j][1]), 1, pF);
    }
  }
  fwrite(&step, sizeof(step), 1, pF);
  fclose(pF);

  /* close "probes.txt" file */
  fclose(pFprobest);

  /* close "probes.bin" file */
  fclose(pFprobesb);

} /* end Finalize() */


/*********
*  MAIN  *   Main function
*********/
int main()
{
  /*--- A handler for a SIGINT for the whole program ---*/
  signal(SIGINT, termHandler);

	/*--- Initialization ---*/
  Initialize();

  /*--- Time stepping ---*/
  while (1) {

    /*--- Check exit conditions ---*/
    if (step == numstep || !continFlag)
      break;

    printf("%d\n", ++step);

    /*=== First stage ===*/
    /*--- BC in ghost cells ---*/
    BounCondInGhostCells();
    /*--- Parameters at the cell boundaries ---*/
    Reconstruction();
    /*--- BC at cell boundaries at the domain boundary ---*/
    BounCondOnInterfaces();
    /*--- Fluxes ---*/
    Fluxes();
    /*--- Evolution ---*/
    Evolution();
    /*--- Probes 1st stage ---*/
    Probes();

    /*=== Second stage ===*/
    /*--- BC in ghost cells (2) ---*/
    BounCondInGhostCells();
    /*--- Parameters at the cell boundaries ---*/
    Reconstruction();
    /*--- BC at cell boundaries at the domain boundary ---*/
    BounCondOnInterfaces();
    /*--- Fluxes ---*/
    Fluxes();
    /*--- Evolution ---*/
    Evolution();
    /*--- Probes 2nd stage ---*/
    Probes();

    /*--- Taking frame ---*/
    if (step % f_step == 0)
    {
        //Output();
        Output_VTK();
    }


    /*--- Debug ---*/
    /*if (step == 1)
        Output_VTK();*/

  } /* end while */

  /*--- Output of the flowfield ---*/
  //Output();
  Output_VTK();

  /*--- Backup solution ---*/
  Finalize();

  fprintf(stdout, "simulation stopped!\n");

  return 0;

} /* end main() */

/*-- end duct.c --*/
