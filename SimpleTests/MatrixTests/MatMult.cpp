#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <omp.h>
#include <cstdlib>
#include "sparse.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/CRSMatrixOpenMP.h"
#include "LinAlg/Sparse/CRSMatrixPThreads.h"
#include "RunTimeTimer.h"
#include "CPUTimeTimer.h"

int main(int argc, char *argv[])
{
	if (argc < 4) {
		std::cout << "Usage: " << argv[0] << " num_of_threads matrix number_of_multiplications resultfile" << std::endl;
		exit (1);
	}

	// read number of threads
	unsigned n_threads (1);
	n_threads = atoi (argv[1]);

	// read the number of multiplication to execute
	unsigned n_mults (0);
	n_mults = atoi (argv[3]);

	std::string fname_mat (argv[2]);

	bool verbose (true);

	// *** reading matrix in crs format from file
	std::ifstream in(fname_mat.c_str(), std::ios::in | std::ios::binary);
	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		if (verbose) {
			std::cout << "reading matrix from " << fname_mat << " ... " << std::flush;
		}
		RunTimeTimer timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		timer.stop();
		if (verbose) {
			std::cout << "ok, " << timer.elapsed() << " s)" << std::endl;
		}
	} else {
		std::cout << "error reading matrix from " << fname_mat << std::endl;
	}
	unsigned nnz(iA[n]);
	if (verbose) {
		std::cout << "Parameters read: n=" << n << ", nnz=" << nnz << std::endl;
	}

/*
	// test: fill all nnz entries with 1
	for (unsigned k(0); k<nnz; k++) A[k] = 2.0;
*/

	MathLib::CRSMatrix<double> mat (n, iA, jA, A);
	std::cout << mat.getNRows() << " x " << mat.getNCols() << std::endl;
//	CRSMatrixOpenMP<double> mat (n, iA, jA, A, n_threads);
//	CRSMatrixPThreads<double> mat (n, iA, jA, A, n_threads);

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	if (verbose) {
		std::cout << "matrix vector multiplication with Toms amuxCRS (" << n_threads << " threads) ... " << std::flush;
	}
	RunTimeTimer run_timer;
	CPUTimeTimer cpu_timer;
	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
//		std::cout << "mult " << k << " ... " << std::flush;
		mat.amux (1.0, x, y);
//		std::cout << "ok" << std::endl;
	}
	cpu_timer.stop();
	run_timer.stop();

/*
	// test result of matrix vector mult
	std::cout << "simple test of matrix vector mult ... " << std::flush;
	for (unsigned k(0); k<n; k++) {
		if (std::fabs(y[k] - 2*(iA[k+1]-iA[k])) > std::numeric_limits<double>::min()) {
			std::cout << "value y[" << k << "] = " << y[k] << " not equal to number of entries in row " << k << ", that is " << iA[k+1]-iA[k] << std::endl;
		}
	}
*/

	if (verbose) {
		std::cout << "done ["
<< cpu_timer.elapsed() << " sec cpu time], [" << run_timer.elapsed() << " sec run time]" << std::endl;
	} else {
		if (argc == 5) {
			std::ofstream result_os (argv[4], std::ios::app);
			if (result_os) {
				result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
			}
			result_os.close();
		} else {
			std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
		}
	}

	delete [] x;
	delete [] y;

	return 0;
}

