/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \qfile MatVecMultNDPermOpenMP.cpp
 *
 * Created on 2012-01-20 by Thomas Fischer
 */

#include <cstdlib>

// BaseLib
#include "RunTime.h"
#include "CPUTime.h"
// BaseLib/tclap
#include "tclap/CmdLine.h"
// ThirdParty/logog
#include "logog/include/logog.hpp"
#include "logog/include/formatter.hpp"

// MathLib
#include "sparse.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReorderedOpenMP.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"

#ifdef UNIX
#include <sys/unistd.h>
#endif

#ifdef OGS_BUILD_INFO
#include "BuildInfo.h"
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * new formatter for logog
 */
class FormatterCustom : public logog::FormatterGCC
{
    virtual TOPIC_FLAGS GetTopicFlags( const logog::Topic &topic )
    {
        return ( Formatter::GetTopicFlags( topic ) &
                 ~( TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG ));
    }
};

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();

	TCLAP::CmdLine cmd("The purpose of this program is the speed test of sparse matrix vector multiplication (MVM) employing OpenMP technique, where the matrix is stored in CRS format. Before executing the MVM a nested dissection reordering is performed.", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m matrix".
	TCLAP::ValueArg<std::string> matrix_arg("m","matrix","input matrix file in CRS format",true,"","file name of the matrix in CRS format");

	// Add the argument matrix_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( matrix_arg );

	TCLAP::ValueArg<unsigned> n_cores_arg("p", "number-cores", "number of cores to use", true, 1, "number of cores");
	cmd.add( n_cores_arg );

	TCLAP::ValueArg<unsigned> n_mults_arg("n", "number-of-multiplications", "number of multiplications to perform", true, 10, "number of multiplications");
	cmd.add( n_mults_arg );

	TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
	cmd.add( output_arg );

	TCLAP::ValueArg<bool> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "string");
	cmd.add( verbosity_arg );

	cmd.parse( argc, argv );

	// read the number of multiplication to execute
	unsigned n_mults (n_mults_arg.getValue());
	std::string fname_mat (matrix_arg.getValue());
	bool verbose (verbosity_arg.getValue());

	FormatterCustom *custom_format (new FormatterCustom);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	// read number of threads
	unsigned n_threads (n_cores_arg.getValue());

#ifdef OGS_BUILD_INFO
	INFO("%s was build with compiler %s", argv[0], CMAKE_CXX_COMPILER);
#ifdef NDEBUG
	INFO("CXX_FLAGS: %s %s", CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_RELEASE);
#else
	INFO("CXX_FLAGS: %s %s", CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_DEBUG);
#endif
#endif

#ifdef UNIX
	const size_t length(256);
	char *hostname(new char[length]);
	gethostname (hostname, length);
	INFO("hostname: %s", hostname);
	delete [] hostname;
#endif

	// *** reading matrix in crs format from file
	std::ifstream in(fname_mat.c_str(), std::ios::in | std::ios::binary);
	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		if (verbose) {
			INFO("reading matrix from %s ...", fname_mat.c_str());
		}
		BaseLib::RunTime timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		timer.stop();
		if (verbose) {
			INFO("\t- took %e s", timer.elapsed());
		}
	} else {
		ERR("error reading matrix from %s", fname_mat.c_str());
		return -1;
	}
	unsigned nnz(iA[n]);
	if (verbose) {
		INFO("\tParameters read: n=%d, nnz=%d", n, nnz);
	}

#ifdef _OPENMP
	omp_set_num_threads(n_threads);
	MathLib::CRSMatrixReorderedOpenMP mat(n, iA, jA, A);
#else
	delete [] iA;
	delete [] jA;
	delete [] A;
	ERROR("program is not using OpenMP");
	return -1;
#endif
	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	// create time measurement objects
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;

	// calculate the nested dissection reordering
	if (verbose) {
		INFO("*** calculating nested dissection (ND) permutation of matrix ...");
	}
	run_timer.start();
	cpu_timer.start();
	MathLib::Cluster cluster_tree(n, iA, jA);
	unsigned *op_perm(new unsigned[n]);
	unsigned *po_perm(new unsigned[n]);
	for (unsigned k(0); k<n; k++)
		op_perm[k] = po_perm[k] = k;
	cluster_tree.createClusterTree(op_perm, po_perm, 1000);
	cpu_timer.stop();
	run_timer.stop();
	if (verbose) {
		INFO("\t[ND] - took %e sec \t%e sec", cpu_timer.elapsed(), run_timer.elapsed());
	}

	// applying the nested dissection reordering
	if (verbose) {
		INFO("\t[ND] applying nested dissection permutation to FEM matrix ... ");
	}
	run_timer.start();
	cpu_timer.start();
	mat.reorderMatrix(op_perm, po_perm);
	cpu_timer.stop();
	run_timer.stop();
	if (verbose) {
		INFO("\t[ND]: - took %e sec\t%e sec", cpu_timer.elapsed(), run_timer.elapsed());
	}

	if (verbose) {
		INFO("*** %d matrix vector multiplications (MVM) with Toms amuxCRS (%d threads)... ", n_mults, n_threads);
	}

	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
		mat.amux (1.0, x, y);
	}
	cpu_timer.stop();
	run_timer.stop();

	if (verbose) {
		INFO("\t[MVM] - took %e sec cpu time, %e sec run time", cpu_timer.elapsed(), run_timer.elapsed());
	}

	delete [] x;
	delete [] y;

	delete custom_format;
	delete logogCout;
	LOGOG_SHUTDOWN();

	return 0;
}
