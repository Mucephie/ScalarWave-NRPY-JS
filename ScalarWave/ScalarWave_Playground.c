// Part P0: Set the number of ghost cells, from NRPy+'s FD_CENTDERIVS_ORDER
#define NGHOSTS 2

const int NSKIP_2D_OUTPUT = 5;

// Part P1: Import needed header files
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include <emscripten.h>

// Part P2: Add needed #define's to set data type, the IDX4() macro, and the gridfunctions
// Part P2a: set REAL=double, so that all floating point numbers are stored to at least ~16 significant digits.
#define REAL double
// Part P2b: Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//           data in a 1D array. In this case, consecutive values of "i" 
//           (all other indices held to a fixed value) are consecutive in memory, where 
//           consecutive values of "j" (fixing all other indices) are separated by 
//           Nxx_plus_2NGHOSTS[0] elements in memory. Similarly, consecutive values of
//           "k" are separated by Nxx_plus_2NGHOSTS[0]*Nxx_plus_2NGHOSTS[1] in memory, etc.
#define IDX4(g,i,j,k) \
( (i) + Nxx_plus_2NGHOSTS[0] * ( (j) + Nxx_plus_2NGHOSTS[1] * ( (k) + Nxx_plus_2NGHOSTS[2] * (g) ) ) )
// Part P2c: Set UUGF and VVGF macros
#define NUM_GFS 2
#define UUGF 0
#define VVGF 1

// Step P3: Set free parameters for the initial data
const REAL wavespeed = 1.0;
const REAL kk0 = 1.0;
const REAL kk1 = 1.0;
const REAL kk2 = 1.0;


// Part P5: Declare the function to evaluate the scalar wave RHSs
void rhs_eval(const int Nxx[3], const int Nxx_plus_2NGHOSTS[3], const REAL dxx[3], const REAL *in_gfs, REAL *rhs_gfs) {
#include "ScalarWave_RHSs.h"
}


// Part P6: Declare boundary condition FACE_UPDATE macro,
//          which updates a single face of the 3D grid cube
//          using quadratic polynomial extrapolation.
const int MAXFACE = -1;
const int NUL = +0;
const int MINFACE = +1;
#define  FACE_UPDATE(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
        gfs[IDX4(which_gf,i0,i1,i2)] =                                  \
          +3.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -3.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]  \
          +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
      }


// Part P7: Boundary condition driver routine: Apply BCs to all six
//          boundary faces of the cube, filling in the innermost
//          ghost zone first, and moving outward.
void apply_bcs(const int Nxx[3], const int Nxx_plus_2NGHOSTS[3], REAL *gfs) {
#pragma omp parallel for
	for (int which_gf = 0; which_gf < NUM_GFS; which_gf++) {
		int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
		int imax[3] = { Nxx_plus_2NGHOSTS[0] - NGHOSTS, Nxx_plus_2NGHOSTS[1] - NGHOSTS, Nxx_plus_2NGHOSTS[2] - NGHOSTS };
		for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
			// After updating each face, adjust imin[] and imax[] 
			//   to reflect the newly-updated face extents.
			FACE_UPDATE(which_gf, imin[0] - 1, imin[0], imin[1], imax[1], imin[2], imax[2], MINFACE, NUL, NUL); imin[0]--;
			FACE_UPDATE(which_gf, imax[0], imax[0] + 1, imin[1], imax[1], imin[2], imax[2], MAXFACE, NUL, NUL); imax[0]++;

			FACE_UPDATE(which_gf, imin[0], imax[0], imin[1] - 1, imin[1], imin[2], imax[2], NUL, MINFACE, NUL); imin[1]--;
			FACE_UPDATE(which_gf, imin[0], imax[0], imax[1], imax[1] + 1, imin[2], imax[2], NUL, MAXFACE, NUL); imax[1]++;

			FACE_UPDATE(which_gf, imin[0], imax[0], imin[1], imax[1], imin[2] - 1, imin[2], NUL, NUL, MINFACE); imin[2]--;
			FACE_UPDATE(which_gf, imin[0], imax[0], imin[1], imax[1], imax[2], imax[2] + 1, NUL, NUL, MAXFACE); imax[2]++;
		}
	}
}


// Get the evolved grid functions
// Should be complete
//
EMSCRIPTEN_KEEPALIVE
REAL * get_gfs() {
	return evol_gfs;
}

// Get size of simulation space 
// Search NUM_EVOL_GFS and check setup of Ntot
// Could have to do with NUM_GFS in step 0c
//
// EMSCRIPTEN_KEEPALIVE
// int get_gfs_size(){
//
//
// return NUM_EVOL_GFS * Ntot;
//}



// Initialize simulation 
// Needs initial data from exact(t=0), or random data like line 187 SENR-em.c
//
EMSCRIPTEN_KEEPALIVE
void init_sim() {
	const int argv = 128;
	// Step 0b: Set up numerical grid structure, first in space...
	const int Nx0x1x2 = atoi(argv); // What does atoi() do?
	const int Nxx[3] = { Nx0x1x2, Nx0x1x2, Nx0x1x2 };
	const int Nxx_plus_2NGHOSTS[3] = { Nxx[0] + 2 * NGHOSTS, Nxx[1] + 2 * NGHOSTS, Nxx[2] + 2 * NGHOSTS };
	const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS[0] * Nxx_plus_2NGHOSTS[1] * Nxx_plus_2NGHOSTS[2];

	const REAL xxmin[3] = { -10.,-10.,-10. };
	const REAL xxmax[3] = { 10., 10., 10. };


	//          ... and then set up the numerical grid structure in time:
	const REAL t_final = xxmax[0] * 0.8; /* Final time is set so that at t=t_final,
										  data at the origin have not been corrupted
										  by the approximate outer boundary condition */
	const REAL CFL_FACTOR = 0.5; // Set the CFL Factor

	// Step 0c: Allocate memory for gridfunctions
	REAL *evol_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
	REAL *next_in_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
	REAL *k1_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
	REAL *k2_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
	REAL *k3_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);
	REAL *k4_gfs = (REAL *)malloc(sizeof(REAL) * NUM_GFS * Nxx_plus_2NGHOSTS_tot);

	// Step 0d: Set up coordinates: Set dx, and then dt based on dx_min and CFL condition
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
// xx[0][i] = xxmin[0] + (i-NGHOSTS)*dxx[0]
	REAL dxx[3];
	for (int i = 0; i < 3; i++) dxx[i] = (xxmax[i] - xxmin[i]) / ((REAL)Nxx[i]);

	
	int Nt = (int)(t_final / dt + 0.5); // The number of points in time.
										//Add 0.5 to account for C rounding down integers.
	


	// Step 0e: Set up Cartesian coordinate grids
	REAL *xx[3];
	for (int i = 0; i < 3; i++) {
		xx[i] = (REAL *)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS[i]);
		for (int j = 0; j < Nxx_plus_2NGHOSTS[i]; j++) {
			xx[i][j] = xxmin[i] + (j - NGHOSTS)*dxx[i];
		}
	}

	// Step 1: Set up initial data to be exact solution at time=0:
	// 
	// We don't have exact_solution() anymore, will change to take input from mouse. 
	//
	exact_solution(Nxx_plus_2NGHOSTS, 0.0, xx, evol_gfs);
	output_2D(0, 0.0, evol_gfs, evol_gfs, Nxx, Nxx_plus_2NGHOSTS, xx);
}



// Run next step of simulation 
// Configure N_iters
//
EMSCRIPTEN_KEEPALIVE
void run_sim(int iter_start, int N_iters)
{
	int iter_start = 0;
	REAL dt = CFL_FACTOR * MIN(dxx[0], MIN(dxx[1], dxx[2])); // CFL condition

	for (int n = iter_start; n <= Nt; n++) { // Main loop to progress forward in time.

		const REAL t = n * dt;

	  // Step 2b: Evolve scalar wave initial data forward in time using Method of Lines with RK4 algorithm,
	  //          applying quadratic extrapolation outer boundary conditions.

	  /***************************************************/
	  /* Implement RK4 for Method of Lines timestepping: */
	  /***************************************************/


	  /* -= RK4: Step 1 of 4 =- */
	  /* First evaluate k1 = RHSs expression             */
		rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, evol_gfs, k1_gfs);
		/* Next k1 -> k1*dt, and then set the input for    */
		/*    the next RHS eval call to y_n+k1/2           */
		for (int i = 0; i < Nxx_plus_2NGHOSTS_tot*NUM_GFS; i++) {
			k1_gfs[i] *= dt;
			next_in_gfs[i] = evol_gfs[i] + k1_gfs[i] * 0.5;
		}
		/* Finally, apply boundary conditions to           */
		/* next_in_gfs, so its data are set everywhere.    */
		apply_bcs(Nxx, Nxx_plus_2NGHOSTS, next_in_gfs);

		/* -= RK4: Step 2 of 4 =- */
		rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, next_in_gfs, k2_gfs);
		for (int i = 0; i < Nxx_plus_2NGHOSTS_tot*NUM_GFS; i++) {
			k2_gfs[i] *= dt;
			next_in_gfs[i] = evol_gfs[i] + k2_gfs[i] * 0.5;
		}
		apply_bcs(Nxx, Nxx_plus_2NGHOSTS, next_in_gfs);

		/* -= RK4: Step 3 of 4 =- */
		rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, next_in_gfs, k3_gfs);
		for (int i = 0; i < Nxx_plus_2NGHOSTS_tot*NUM_GFS; i++) {
			k3_gfs[i] *= dt;
			next_in_gfs[i] = evol_gfs[i] + k3_gfs[i];
		}
		apply_bcs(Nxx, Nxx_plus_2NGHOSTS, next_in_gfs);

		/* -= RK4: Step 4 of 4 =- */
		rhs_eval(Nxx, Nxx_plus_2NGHOSTS, dxx, next_in_gfs, k4_gfs);
		for (int i = 0; i < Nxx_plus_2NGHOSTS_tot*NUM_GFS; i++) {
			k4_gfs[i] *= dt;

			evol_gfs[i] += (1.0 / 6.0)*(k1_gfs[i] + 2.0*k2_gfs[i] + 2.0*k3_gfs[i] + k4_gfs[i]);
		}
		apply_bcs(Nxx, Nxx_plus_2NGHOSTS, evol_gfs);
	}

}

// Clean things up/ clear memory of this step
// Should be complete now
//
EMSCRIPTEN_KEEPALIVE
void clean_sim() {
	// Step : Free all allocated memory
	free(k4_gfs);
	free(k3_gfs);
	free(k2_gfs);
	free(k1_gfs);
	free(next_in_gfs);
	free(evol_gfs);
	for (int i = 0; i < 3; i++) free(xx[i]);
}