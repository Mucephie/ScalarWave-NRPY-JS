/**
 * Copy of SENR.c, but with file I/O removed,
 * and run arguments explicitly set.
 */
#include "SENR.h"
#include <emscripten.h>

EMSCRIPTEN_KEEPALIVE
void print_logo()
{
	// Print the SENR logo.
#include "Logo.c"
}

// Global variables
paramstruct params;
REAL *x1G, *x2G, *x3G;
REAL *gfs_n, *gfs_1, *gfs_np1, *gfs_k, *gfs_aux, *yy;
precomputed_quantities *precomp;
// Physical coordinate bounds
REAL y_min[DIM], y_max[DIM];


EMSCRIPTEN_KEEPALIVE
REAL * get_gfs()
{
	return gfs_n;
}

EMSCRIPTEN_KEEPALIVE
int get_gfs_size()
{
	const int Nx1 = params.Nx1;
	const int Nx2 = params.Nx2;
	const int Nx3 = params.Nx3;
	const int NGHOSTS = params.NGHOSTS;
	params.Npts1 = (Nx1 + 2 * NGHOSTS);
	params.Npts2 = (Nx2 + 2 * NGHOSTS);
	params.Npts3 = (Nx3 + 2 * NGHOSTS);
	const int Npts1 = params.Npts1;
	const int Npts2 = params.Npts2;
	const int Npts3 = params.Npts3;
	const int Ntot = Npts1 * Npts2 * Npts3;
	return NUM_EVOL_GFS * Ntot;
}


EMSCRIPTEN_KEEPALIVE
void init_sim()
{

	// emscripten: explicitly provide arguments here
	int argc = 6;
	const char* argv[] = { "SENR", "128", "32", "2", "1000.0", "0.15" };

	// Error out if command-line arguments are not set. Otherwise set these "runtime parameters"
#include "runtime_params_set-NRPyGEN.h"

  // Set basic simulation parameters, including setting floating point precision
  // (double or long double), coord system, ID type, etc.
#include "Input_CompileTime_params-NRPyGEN.h"

  // Display simulation info
	printf("\x1B[32mID: %s, Evol: %s, Coords: %s, FD order: %d\x1B[0m\n",
		params.ID_scheme, params.Evol_scheme, params.CoordSystem, params.FDCENTERDERIVS_FDORDER);

#ifdef USE_LONG_DOUBLE
	if (strncmp(params.PRECISION, "long double", 100) != 0) {
		printf("ERROR: You set USE_LONG_DOUBLE in loop_setup.h, but forgot to tell NRPy!\n");
		exit(1);
	}
#else
	if (strncmp(params.PRECISION, "long double", 100) == 0) {
		printf("ERROR: You told NRPy to use long doubles, but forgot to set USE_LONG_DOUBLE in loop_setup.h!\n");
		exit(1);
	}
#endif

	const int FDCENTERDERIVS_FDORDER = params.FDCENTERDERIVS_FDORDER;

	params.Nx1 = atoi(argv[1]); params.Nx2 = atoi(argv[2]); params.Nx3 = atoi(argv[3]);
	params.NGHOSTS = ((params.FDCENTERDERIVS_FDORDER) / 2 + 1);

	if (params.UPWIND == 0 || params.UPWIND == -1) {
		params.NGHOSTS = ((params.FDCENTERDERIVS_FDORDER) / 2);
	}

	const int Nx1 = params.Nx1;
	const int Nx2 = params.Nx2;
	const int Nx3 = params.Nx3;
	params.del[0] = (params.x1max - params.x1min) / ((REAL)Nx1);
	params.del[1] = (params.x2max - params.x2min) / ((REAL)Nx2);
	params.del[2] = (params.x3max - params.x3min) / ((REAL)Nx3);

	const int NGHOSTS = params.NGHOSTS;
	params.Npts1 = (Nx1 + 2 * NGHOSTS);
	params.Npts2 = (Nx2 + 2 * NGHOSTS);
	params.Npts3 = (Nx3 + 2 * NGHOSTS);

	const int Npts1 = params.Npts1;
	const int Npts2 = params.Npts2;
	const int Npts3 = params.Npts3;

	// Total number of grid points, including ghost zones
	const int Ntot = Npts1 * Npts2 * Npts3;

	const double del[3] = { params.del[0],params.del[1],params.del[2] };
	const double x123min[3] = { params.x1min,params.x2min,params.x3min };

	// Uniform coordinates
#define MIN(i,j) ((i < j) ? i : j)
#define MAX(i,j) ((i > j) ? i : j)
	x1G = (REAL *)malloc(sizeof(REAL)*Npts1);
	x2G = (REAL *)malloc(sizeof(REAL)*Npts2);
	x3G = (REAL *)malloc(sizeof(REAL)*Npts3);

	// Initialize static (uniform) cell-centered coordinates on the chosen coordinate domain.
	//            Add 0.5 to enforce cell-centering.
	for (int ii = 0; ii < Npts1; ii++) x1G[ii] = x123min[0] + del[0] * ((REAL)(ii - NGHOSTS) + 0.5);
	for (int jj = 0; jj < Npts2; jj++) x2G[jj] = x123min[1] + del[1] * ((REAL)(jj - NGHOSTS) + 0.5);
	for (int kk = 0; kk < Npts3; kk++) x3G[kk] = x123min[2] + del[2] * ((REAL)(kk - NGHOSTS) + 0.5);

	// Static coordinate domains set by the grid resolution
	const REAL x_min[3] = { x1G[NGHOSTS], x2G[NGHOSTS], x3G[NGHOSTS] };
	const REAL x_max[3] = { x1G[Npts1 - NGHOSTS - 1], x2G[Npts2 - NGHOSTS - 1], x3G[Npts3 - NGHOSTS - 1] };

#ifdef AUTO_BOUNDARY
	// SymmetryMap and ParityMap for assigning ghost zone values
	int *SymmetryMap = (int *)malloc(sizeof(int)*Ntot);
	REAL *ParityMap = (REAL *)malloc(sizeof(REAL)*NUM_EVOL_GFS*Ntot);

	// Find the ghost zone map
	Map_Ghosts(x1G, x2G, x3G, SymmetryMap, ParityMap, params);
#endif

	// Allocate grid function storage

	// Evolved fields
	gfs_n = (REAL *)malloc(sizeof(REAL)*NUM_EVOL_GFS*Ntot);
	gfs_1 = (REAL *)malloc(sizeof(REAL)*NUM_EVOL_GFS*Ntot);
	gfs_np1 = (REAL *)malloc(sizeof(REAL)*NUM_EVOL_GFS*Ntot);

	// RK4 operators
	gfs_k = (REAL *)malloc(sizeof(REAL)*NUM_EVOL_GFS*Ntot);

	// Auxiliary grid functions
	gfs_aux = (REAL *)malloc(sizeof(REAL)*NUM_AUX_GFS*Ntot);

	// Curvilinear coordinates
	yy = (REAL *)malloc(sizeof(REAL)*DIM*Ntot); // (y1, y2, y3)

	// Poison the grid functions with NaNs
#ifdef POISON_GFS
	LOOP_GZFILL(ii, jj, kk) {
		const int idx = IDX3(ii, jj, kk);
		for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
			gfs_n[IDX4pt(gf, idx)] = 1.0 / 0.0;
			gfs_1[IDX4pt(gf, idx)] = 1.0 / 0.0;
			gfs_np1[IDX4pt(gf, idx)] = 1.0 / 0.0;
			gfs_k[IDX4pt(gf, idx)] = 1.0 / 0.0;
		}
		// AUX gridfunctions
		for (int gf = 0; gf < NUM_AUX_GFS; gf++) {
			gfs_aux[IDX4pt(gf, idx)] = 1.0 / 0.0;
		}
	}
#endif

	precomp = (precomputed_quantities *)malloc(sizeof(precomputed_quantities)*Ntot);
	// Store coordinate range
	LOOP_GZFILL(ii, jj, kk) {
		const int idx = IDX3(ii, jj, kk);

		const REAL x1 = x1G[ii];
		const REAL x2 = x2G[jj];
		const REAL x3 = x3G[kk];

		{
#include "reference_metric/NRPy_codegen/NRPy_compute__precomputed_hatted_quantities.h" // Located in ../common_functions/
		}

		{
#include "reference_metric/NRPy_codegen/yy.h" // Located in ../common_functions/
		}
	}

#ifdef ENABLE_RANDOM_ID
	// Fill grid functions with random initial data
	Set_Random_Initial_Data(gfs_n, params);
#else
	// Set BSSN initial data from NRPy
	Set_Initial_Data(x1G, x2G, x3G, gfs_n, params);
#endif

	// Set up spatially-varying Kreiss-Oliger dissipation strength gridfunction
	LOOP_GZFILL(i, j, k) {
		const int idx = IDX3(i, j, k);
		const REAL x1 = x1G[i];
		const REAL x2 = x2G[j];
		const REAL x3 = x3G[k];
#include "reference_metric/NRPy_codegen/dist_from_origin.h" // Located in ../common_functions/
		gfs_aux[IDX4pt(KOSTRNGH, idx)] = KreissOligerDissipation_vs_DistFromOrigin(dist_from_origin);
	}

	// Apply boundary conditions on the initial data or data read from checkpoint.
#pragma omp parallel for
	for (int gf = 0; gf < NUM_EVOL_GFS; gf++) {
		Apply_out_sym_bcs(gf, 0, x1G, x2G, x3G, gfs_n, params);
	}

	// Evaluate detg on the initial data or data read from checkpoint.
	BSSN_Det(x1G, x2G, x3G, gfs_n, gfs_aux, params);

}


EMSCRIPTEN_KEEPALIVE
void run_sim(int N_iters)
{

	int iter_start = 0;
	const int NGHOSTS = params.NGHOSTS;
	const int Npts1 = params.Npts1;
	const int Npts2 = params.Npts2;
	const int Npts3 = params.Npts3;
	const int Nx1 = params.Nx1;
	const int Nx2 = params.Nx2;
	const int Nx3 = params.Nx3;
	const REAL x_min[3] = { x1G[NGHOSTS], x2G[NGHOSTS], x3G[NGHOSTS] };
	const REAL x_max[3] = { x1G[Npts1 - NGHOSTS - 1], x2G[Npts2 - NGHOSTS - 1], x3G[Npts3 - NGHOSTS - 1] };
	const double del[3] = { params.del[0],params.del[1],params.del[2] };

	// Start the timer, for benchmarking
	struct timespec start, end;
	clock_gettime(CLOCK_REALTIME, &start);
	const int timestamp = (int)start.tv_sec;


	// Time step, satisfying the CFL condition
	dt = CFL * Find_Time_Step(x1G, x2G, x3G, yy, gfs_n, gfs_aux, y_min, y_max, params);

	// Output frame index
	int frame_ind = 0;

	// Display grid info
	for (int dirn = 0; dirn < DIM; dirn++) {
		printf("x%dmin = %f, x%dmax = %f, Dx%d = %.2e\n", dirn, (double)x_min[dirn], dirn,
			(double)x_max[dirn], dirn, (double)del[dirn]);
		printf("y%dmin = %f, y%dmax = %f\n", dirn, (double)y_min[dirn], dirn, (double)y_max[dirn]);
	}

	// Iterate through time

	for (int n = iter_start; n <= N_iters; n++) {
		// Current time
		const REAL t = n * dt;

		////////////////////////////////////////////////////////////////////////
		// Numerical integration
#include "evolution_equations/RK4_Steps.c"

	// Rescale gammabarDD to "unit" determinant
		Rescale_Metric_Det(x1G, x2G, x3G, gfs_np1, gfs_aux, params);

		// Enforce trace-free condition on AbarDD
		//Remove_Trace(x1G,x2G,x3G, yy, gfs_np1, params);

		////////////////////////////////////////////////////////////////////////

		int exit_flag = 0;

		// Update time step, by simply swapping pointers in memory:
		// First declare a temporary pointer to gfs_n:
		REAL *tmp_pointer = gfs_n;
		// u^{n} -> u^{n+1}
		// Then set gfs_n to point to gfs_np1
		gfs_n = gfs_np1;
		// Finally, set gfs_np1 to the temporary pointer.
		gfs_np1 = tmp_pointer;

		// Check for NaNs in the updated gridfunctions
		NaN_Checker(0, n, t, x1G, x2G, x3G, gfs_n, gf_name, params);

		////////////////////////////////////////////////////////////////////////

		// Measure average time per iteration
		clock_gettime(CLOCK_REALTIME, &end);
		const long long unsigned int time_in_ns = 1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
		const REAL s_per_iteration_avg = ((REAL)time_in_ns / (REAL)n) / 1.0e9;

		const int iterations_remaining = N_iters - n;
		const REAL time_remaining_in_mins = s_per_iteration_avg * (REAL)iterations_remaining / 60.0;

		const REAL num_RHS_pt_evals = (REAL)(Nx1*Nx2*Nx3) * 4.0 * (REAL)n; // 4 RHS evals per gridpoint for RK4
		const REAL RHS_pt_evals_per_sec = num_RHS_pt_evals / ((REAL)time_in_ns / 1.0e9);

		// Progress indicator printing to stdout
		printf("%c[2K", 27); // Clear the line
		printf("It: %d t=%.2f | %.1f%%; ETA %.0f m | t/h %.2f | gp/s %.2e",
			n, n * (double)dt, (double)(100.0 * (REAL)n / (REAL)N_iters),
			(double)time_remaining_in_mins, (double)(dt * 3600.0 / s_per_iteration_avg), (double)RHS_pt_evals_per_sec);
		fflush(stdout); // Flush the stdout buffer
	} // END for(n <= N_iters)

	printf("\n");
}


void clean_sim()
{
	// Free allocated memory
	free(gfs_n);
	free(gfs_1);
	free(gfs_np1);

	// RK4 operators
	free(gfs_k);

	free(gfs_aux);

	free(yy);
	free(x1G);
	free(x2G);
	free(x3G);

#ifdef AUTO_BOUNDARY
	free(SymmetryMap);
	free(ParityMap);
#endif

	free(precomp);
}
