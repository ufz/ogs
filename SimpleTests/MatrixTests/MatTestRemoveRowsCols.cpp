/**
 * \file
 * \author Thomas Fischer
 * \date   2011-11-08
 * \brief  Test for removing entries in matrices.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <fstream>
#include <iostream>

// BaseLib
#include "RunTime.h"
#include "CPUTime.h"

// MathLib
#include "LinAlg/Sparse/CRSMatrix.h"

int main(int argc, char *argv[])
{
	if (argc < 3) {
		std::cout << "Usage: " << argv[0] << " input-matrix output-matrix" << std::endl;
		return 1;
	}

	std::string fname_mat (argv[1]);
	bool verbose (true);

	// *** reading matrix in crs format from file
	std::ifstream in(fname_mat.c_str(), std::ios::in | std::ios::binary);
	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		if (verbose) {
			std::cout << "reading matrix from " << fname_mat << " ... " << std::flush;
		}
		BaseLib::RunTime timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		in.close();
		timer.stop();
		if (verbose) {
			std::cout << "ok, " << timer.elapsed() << " s" << std::endl;
		}
	} else {
		std::cout << "error reading matrix from " << fname_mat << std::endl;
	}
	unsigned nnz(iA[n]);
	if (verbose) {
		std::cout << "Parameters read: n=" << n << ", nnz=" << nnz << std::endl;
	}

	MathLib::CRSMatrix<double, unsigned> *mat (new MathLib::CRSMatrix<double, unsigned>(n, iA, jA, A));

	const unsigned n_rows_cols_to_erase(300);
	unsigned *rows_cols_to_erase(new unsigned[n_rows_cols_to_erase]);

	for (unsigned k(0); k<n_rows_cols_to_erase; k++) {
		rows_cols_to_erase[k] = (k+1)*2;
	}

	BaseLib::RunTime timer;
	std::cout << "erasing " << n_rows_cols_to_erase << " rows and columns ... " << std::flush;
	timer.start();
	mat->eraseEntries(n_rows_cols_to_erase, rows_cols_to_erase);
	timer.stop();
	std::cout << "ok, " << timer.elapsed() << " s" << std::endl;
	delete[] rows_cols_to_erase;

	fname_mat = argv[2];
	std::ofstream out (fname_mat.c_str(), std::ios::binary);
	CS_write (out, mat->getNRows(), mat->getRowPtrArray(), mat->getColIdxArray(), mat->getEntryArray());
	out.close();

	std::cout << "wrote " << fname_mat << " with " << mat->getNRows() << " rows and " << mat->getRowPtrArray()[mat->getNRows()] << " entries" << std::endl;

	delete mat;
}
