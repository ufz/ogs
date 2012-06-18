#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include "sparse.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/CRSMatrixOpenMP.h"
#include "LinAlg/Sparse/CRSMatrixPThreads.h"
#include "RunTime.h"
#include "CPUTime.h"
#include "logog.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;

	if (argc < 4) {
		std::cout << "Usage: " << argv[0] << " num_of_threads matrix number_of_multiplications resultfile" << std::endl;
		INFO("Usage: %s num_of_threads matrix number_of_multiplications resultfile", argv[0]);
		exit (1);
	}

	// read number of threads
	unsigned n_threads (1);
	n_threads = atoi (argv[1]);

	// read the number of multiplication to execute
	unsigned n_mults (0);
	n_mults = atoi (argv[3]);

	std::string fname_mat (argv[2]);

	// *** reading matrix in crs format from file
	std::ifstream in(fname_mat.c_str(), std::ios::in | std::ios::binary);
	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		DBUG("reading matrix from %s ...", fname_mat.c_str());
		BaseLib::RunTime timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		timer.stop();
		DBUG("ok, %e s", timer.elapsed());
	} else {
		ERR("error reading matrix from %s", fname_mat.c_str());
	}
	unsigned nnz(iA[n]);
	INFO("Parameters read: n=%i, nnz=%i", n, nnz);

#ifdef _OPENMP
	omp_set_num_threads(n_threads);
	MathLib::CRSMatrixOpenMP<double, unsigned> mat (n, iA, jA, A);
#else
	MathLib::CRSMatrix<double, unsigned> mat (n, iA, jA, A);
#endif
//	CRSMatrixPThreads<double> mat (n, iA, jA, A, n_threads);
	INFO("%i x %i", mat.getNRows(), mat.getNCols());

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	DBUG("matrix vector multiplication with Toms amuxCRS (%i threads) ...", n_threads);
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;
	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
		mat.amux (1.0, x, y);
	}
	cpu_timer.stop();
	run_timer.stop();

	DBUG("done [%e sec cpu time], [%e sec run time]", cpu_timer.elapsed(), run_timer.elapsed());
	DBUG("CPU time: %e", cpu_timer.elapsed());
	DBUG("wclock time: %e", run_timer.elapsed());

	if (argc == 5) {
		std::ofstream result_os (argv[4], std::ios::app);
		if (result_os) {
			result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
		}
		result_os.close();
	} else {
		INFO("%e \t %e", cpu_timer.elapsed(), run_timer.elapsed());
	}


	delete [] x;
	delete [] y;

	delete logogCout;
	LOGOG_SHUTDOWN();

	return 0;
}

