#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include "LinAlg/Solvers/BiCGStab.h"
#include "LinAlg/Sparse/CRSMatrixDiagPrecond.h"
#include "sparse.h"
#include "vector_io.h"
#include "RunTimeTimer.h"
#include "CPUTimeTimer.h"

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
		x[k] = 1.0;
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
		std::cout << "solving system with BiCGStab method (diagonal preconditioner) ... " << std::flush;

	double eps (1.0e-6);
	unsigned steps (4000);
	RunTimeTimer run_timer;
	CPUTimeTimer cpu_timer;
	run_timer.start();
	cpu_timer.start();

	MathLib::BiCGStab ((*mat), b, x, eps, steps);

	cpu_timer.stop();
	run_timer.stop();

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

