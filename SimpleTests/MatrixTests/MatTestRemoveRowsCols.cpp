/*
 * MatTestRemoveRowsCols.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: TF
 */

#include <fstream>
#include <iostream>

// Base
#include "RunTimeTimer.h"
#include "CPUTimeTimer.h"

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
		RunTimeTimer timer;
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

	const double val0 ((*mat)(1,1));
	std::cout << val0 << std::endl;

//	double val ((*(const_cast<MathLib::CRSMatrix<double, unsigned> const*>(mat)))(1,1));
//	std::cout << val << std::endl;

	MathLib::CRSMatrix<double, unsigned> const& ref_mat(*mat);
	double val1 (ref_mat(1,1));
	std::cout << val1 << std::endl;

	const unsigned n_rows_cols_to_erase(3);
	unsigned *rows_to_erase(new unsigned[n_rows_cols_to_erase]);
	unsigned *cols_to_erase(new unsigned[n_rows_cols_to_erase]);

	for (unsigned k(0); k<n_rows_cols_to_erase; k++) {
		rows_to_erase[k] = (k+1)*2;
		cols_to_erase[k] = (k+1)*2;
	}

	mat->eraseEntries(n_rows_cols_to_erase, rows_to_erase, cols_to_erase);
	delete[] rows_to_erase;
	delete[] cols_to_erase;

	fname_mat = argv[2];
	std::ofstream out (fname_mat.c_str(), std::ios::binary);
	CS_write (out, mat->getNRows(), mat->getRowPtrArray(), mat->getColIdxArray(), mat->getEntryArray());
	out.close();

	delete mat;
}
