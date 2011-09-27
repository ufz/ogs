#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>
#include "CG.h"
#include "CRSMatrixDiagPrecond.h"
#include "sparse.h"
#include "vector_io.h"
#include "timeMeasurement.h"

int main(int argc, char *argv[])
{
	if (argc != 4) {
		std::cout << "Usage: " << argv[0] << " matrix rhs number-of-threads" << std::endl;
		return -1;
	}

	// read number of threads       
	unsigned num_omp_threads (1);
	num_omp_threads = atoi (argv[3]);

	// *** reading matrix in crs format from file
	std::string fname(argv[1]);
	CRSMatrixDiagPrecond *mat (new CRSMatrixDiagPrecond(fname, num_omp_threads));

	unsigned n (mat->getNRows());
	bool verbose (true);
	if (verbose)
		std::cout << "Parameters read: n=" << n << std::endl;

	double *x(new double[n]);
	double *b(new double[n]);

	// *** init start vector x
	for (size_t k(0); k<n; k++) {
		x[k] = 0.0;
	}
	// *** read rhs
	fname = argv[2];
	std::ifstream in(fname.c_str());
	if (in) {
		read (in, n, b);
		in.close();
	} else {
		std::cout << "problem reading rhs - initializing b with 1.0" << std::endl;
		for (size_t k(0); k<n; k++) {
			b[k] = 1.0;
		}
	}

	
	if (verbose)
		std::cout << "solving system with PCG method (diagonal preconditioner) ... " << std::flush;
	time_t start_time, end_time;
	time(&start_time);
	double cg_time (cputime(0.0));
	double eps (1.0e-6);
	unsigned steps (4000);
	CG (mat, b, x, eps, steps, num_omp_threads);
	cg_time = cputime(cg_time);
	time(&end_time);
	if (verbose) {
		std::cout << " in " << steps << " iterations" << std::endl;
		std::cout << "\t(residuum is " << eps << ") took " << cg_time << " sec time and " << (end_time-start_time) << " sec" << std::endl;
	} else {
		std::cout << end_time-start_time << std::endl;
	}

	delete mat;
	delete [] x;
	delete [] b;

	return 0;
}

