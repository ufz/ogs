/**
 * \file
 * \author Thomas Fischer
 * \date   2011-10-05
 * \brief  Test for GMRes preconditioner.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>
#include <cstdlib>

#include "BaseLib/CPUTime.h"
#include "BaseLib/RunTime.h"
#include "MathLib/LinAlg/Solvers/GMRes.h"
#include "MathLib/LinAlg/Sparse/CRSMatrixDiagPrecond.h"
#include "MathLib/vector_io.h"

int main(int argc, char *argv[])
{
	if (argc != 3) {
		std::cout << "Usage: " << argv[0] << " matrix rhs" << std::endl;
		return -1;
	}

	// *** reading matrix in crs format from file
	std::string fname(argv[1]);
	MathLib::CRSMatrixDiagPrecond *mat(new MathLib::CRSMatrixDiagPrecond(fname));

	unsigned n(mat->getNRows());
	bool verbose(true);
	if (verbose) std::cout << "Parameters read: n=" << n << std::endl;

	double *x(new double[n]);
	double *b(new double[n]);

	// *** init start vector x
	for (std::size_t k(0); k < n; k++) {
		x[k] = 1.0;
	}
	// *** read rhs
	fname = argv[2];
	std::ifstream in(fname.c_str());
	if (in) {
		read(in, n, b);
		in.close();
	} else {
		std::cout << "problem reading rhs - initializing b with 1.0"
				<< std::endl;
		for (std::size_t k(0); k < n; k++) {
			b[k] = 1.0;
		}
	}

	if (verbose)
		std::cout << "solving system with GMRes method (diagonal preconditioner) ... "
			<< std::flush;

	double eps(1.0e-6);
	unsigned steps(4000);
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;
	run_timer.start();
	cpu_timer.start();

	MathLib::GMRes((*mat), b, x, eps, 30, steps);

	if (verbose) {
		std::cout << " in " << steps << " iterations" << std::endl;
		std::cout << "\t(residuum is " << eps << ") took "
				<< cpu_timer.elapsed() << " sec time and "
				<< run_timer.elapsed() << " sec" << std::endl;
	} else {
		std::cout << cpu_timer.elapsed() << std::endl;
	}

	delete mat;
	delete[] x;
	delete[] b;

	return 0;
}

