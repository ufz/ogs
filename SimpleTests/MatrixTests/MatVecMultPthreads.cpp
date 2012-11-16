/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MatVecMultPthreads.cpp
 *
 *  Created on  Jul 3, 2012 by Thomas Fischer
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
// ThirdParty/logog
#include "logog/include/logog.hpp"
#include "logog/include/formatter.hpp"
// BaseLib/tclap
#include "tclap/CmdLine.h"

#ifdef UNIX
#include <sys/unistd.h>
#endif

#ifdef OGS_BUILD_INFO
#include "BuildInfo.h"
#endif

#ifdef HAVE_PTHREADS
#include <pthread.h>
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

	TCLAP::CmdLine cmd("Simple matrix vector multiplication test employing pthreads", ' ', "0.1");

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

	TCLAP::ValueArg<bool> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "string");
	cmd.add( verbosity_arg );

	cmd.parse( argc, argv );

	std::string fname_mat (matrix_arg.getValue());

	FormatterCustom *custom_format (new FormatterCustom);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	logog::LogFile *logog_file(NULL);
	if (! output_arg.getValue().empty()) {
		logog_file = new logog::LogFile(output_arg.getValue().c_str());
		logog_file->SetFormatter( *custom_format );
	}

#ifdef OGS_BUILD_INFO
	INFO("%s was build with compiler %s", argv[0], CMAKE_CXX_COMPILER);
	if (std::string(CMAKE_BUILD_TYPE).compare("Release") == 0) {
		INFO("CXX_FLAGS: %s %s", CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_RELEASE);
	} else {
		INFO("CXX_FLAGS: %s %s", CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_DEBUG);
	}
#endif

#ifdef UNIX
	const int max_host_name_len (255);
	char *hostname(new char[max_host_name_len]);
	if (gethostname(hostname, max_host_name_len) == 0)
		INFO("hostname: %s", hostname);
	delete [] host_name_len;
#endif

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
		INFO("\t- took %e s", timer.elapsed());
	} else {
		ERR("error reading matrix from %s", fname_mat.c_str());
		return -1;
	}
	unsigned nnz(iA[n]);
	INFO("\tParameters read: n=%d, nnz=%d", n, nnz);

#ifdef HAVE_PTHREADS
	unsigned n_threads(n_cores_arg.getValue());
	MathLib::CRSMatrixPThreads<double> mat (n, iA, jA, A, n_threads);

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	// read the number of multiplication to execute
	unsigned n_mults (n_mults_arg.getValue());

	INFO("*** %d matrix vector multiplications (MVM) with Toms amuxCRS (%d threads) ...", n_mults, n_threads);
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;
	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
		mat.amux (1.0, x, y);
	}
	cpu_timer.stop();
	run_timer.stop();

	INFO("\t[MVM] - took %e sec cpu time, %e sec run time", cpu_timer.elapsed(), run_timer.elapsed());

	delete [] x;
	delete [] y;
#endif

	delete custom_format;
	delete logogCout;
	delete logog_file;
	LOGOG_SHUTDOWN();

	return 0;
}
