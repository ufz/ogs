/*
 * MatMult.cpp
 *
 *  Created on: Jan 3, 2012
 *      Author: TF
 */

#include <fstream>
#include <iostream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include "sparse.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "LinAlg/Sparse/CRSMatrixOpenMP.h"
#include "LinAlg/Sparse/CRSMatrixPThreads.h"

// BaseLib
#include "RunTime.h"
#include "CPUTime.h"
#include "logog.hpp"
#include "tclap/CmdLine.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();

	TCLAP::CmdLine cmd("Simple matrix vector multiplication test", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m matrix".
	TCLAP::ValueArg<std::string> matrix_arg("m", "matrix", "input matrix file", true, "", "string");

	// Add the argument mesh_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( matrix_arg );

	TCLAP::ValueArg<unsigned> n_cores_arg("p", "number-cores", "number of cores to use", false, 1, "number");
	cmd.add( n_cores_arg );

	TCLAP::ValueArg<unsigned> n_mults_arg("n", "number-of-multiplications", "number of multiplications to perform", true, 10, "number");
	cmd.add( n_mults_arg );

	TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
	cmd.add( output_arg );

	TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "string");
	cmd.add( verbosity_arg );

	cmd.parse( argc, argv );

	// read the number of multiplication to execute
	unsigned n_mults (n_mults_arg.getValue());
	std::string fname_mat (matrix_arg.getValue());

	logog::Cout* logogCout = new logog::Cout;

	// read number of threads
	unsigned n_threads (n_cores_arg.getValue());

	// *** reading matrix in crs format from file
	std::ifstream in(fname_mat.c_str(), std::ios::in | std::ios::binary);
	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		INFO("reading matrix from %s ...", fname_mat.c_str());
		BaseLib::RunTime timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		timer.stop();
		INFO("ok, %e s", timer.elapsed());
	} else {
		ERR("error reading matrix from %s", fname_mat.c_str());
	}
	unsigned nnz(iA[n]);
	INFO("Parameters read: n=%d, nnz=%d", n, nnz);

#ifdef _OPENMP
	omp_set_num_threads(n_threads);
	MathLib::CRSMatrixOpenMP<double, unsigned> mat (n, iA, jA, A);
#else
	MathLib::CRSMatrix<double, unsigned> mat (n, iA, jA, A);
#endif
//	CRSMatrixPThreads<double> mat (n, iA, jA, A, n_threads);
	INFO("%d x %d", mat.getNRows(), mat.getNCols());

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	INFO("matrix vector multiplication with Toms amuxCRS (%d threads) ...", n_threads);
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;
	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
		mat.amux (1.0, x, y);
	}
	cpu_timer.stop();
	run_timer.stop();

	INFO("done [%e sec cpu time], [%e sec run time]", cpu_timer.elapsed(), run_timer.elapsed());
	INFO("CPU time: %e", cpu_timer.elapsed());
	INFO("wclock time: %e", run_timer.elapsed());

	if (! output_arg.getValue().empty()) {
		std::ofstream result_os (output_arg.getValue().c_str(), std::ios::app);
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

