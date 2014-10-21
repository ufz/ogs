#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "LinAlg/Solvers/CG.h"
#include "LinAlg/Sparse/CRSMatrixDiagPrecond.h"
#include "sparse.h"
#include "vector_io.h"
#include "RunTime.h"
#include "CPUTime.h"

#ifdef _OPENMP
#include <omp.h>
#endif

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
	MathLib::CRSMatrixDiagPrecond *mat (new MathLib::CRSMatrixDiagPrecond(fname));

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

	double eps (1.0e-6);
	unsigned steps (4000);
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;
	run_timer.start();
	cpu_timer.start();

	if (num_omp_threads == 1) {
		MathLib::CG(mat, b, x, eps, steps);
	} else {
		#ifdef _OPENMP
		MathLib::CGParallel(mat, b, x, eps, steps);
		#else
		std::cout << "OpenMP is not switched on" << std::endl;
		#endif
	}

	if (verbose) {
		std::cout << " in " << steps << " iterations" << std::endl;
		std::cout << "\t(residuum is " << eps << ") took " << cpu_timer.elapsed() << " sec time and " << run_timer.elapsed() << " sec" << std::endl;
	} else {
		std::cout << cpu_timer.elapsed() << std::endl;
	}

	delete mat;
	delete [] x;
	delete [] b;

	return 0;
}

