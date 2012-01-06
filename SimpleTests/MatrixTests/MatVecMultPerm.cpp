/*
 * MatVecMultPerm.cpp
 *
 *  Created on: Jan 3, 2012
 *      Author: TF
 */

#include <cstdlib>

// BaseLib
#include "RunTimeTimer.h"
#include "CPUTimeTimer.h"

// MathLib
#include "sparse.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReordered.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"
#include "LinAlg/Sparse/CRSMatrix.h"

int main(int argc, char *argv[])
{
	if (argc < 4) {
		std::cout << "Usage: " << argv[0] << " matrix number_of_multiplications resultfile" << std::endl;
		return 1;
	}

	// read the number of multiplication to execute
	unsigned n_mults (0);
	n_mults = atoi (argv[2]);

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

//	MathLib::CRSMatrix<double, unsigned> mat(n, iA, jA, A);
	MathLib::CRSMatrixReordered mat(n, iA, jA, A);
	std::cout << mat.getNRows() << " x " << mat.getNCols() << std::endl;

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	// create time measurement objects
	RunTimeTimer run_timer;
	CPUTimeTimer cpu_timer;

	// calculate the nested dissection reordering
	if (verbose) {
		std::cout << "calculating nested dissection permutation of matrix ... " << std::flush;
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
	if (verbose)
		std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;

	// applying the nested dissection reordering
	if (verbose) {
		std::cout << "applying nested dissection permutation to FEM matrix ... " << std::flush;
	}
	run_timer.start();
	cpu_timer.start();
	mat.reorderMatrix(op_perm, po_perm);
	cpu_timer.stop();
	run_timer.stop();
	if (verbose) std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;

#ifndef NDEBUG
//	MathLib::AdjMat *global_reordered_adj_mat((cluster_tree.getGlobalAdjMat())->getMat(0,n,op_perm, po_perm));
//	const unsigned adj_nnz(global_reordered_adj_mat->getNNZ());
//	double* adj_mat_data(new double[adj_nnz]);
//	for (unsigned k(0); k<adj_nnz; k++) adj_mat_data[k] = 1.0;
//	std::string fname_out (fname_mat);
//	fname_out = fname_out.substr(0,fname_mat.length()-4);
//	fname_out += "_adj.bin";
//	std::ofstream os (fname_out.c_str(), std::ios::binary);
//	CS_write(os, n, global_reordered_adj_mat->getRowPtrArray(), global_reordered_adj_mat->getColIdxArray(), adj_mat_data);
	std::string fname_fem_out (fname_mat);
	fname_fem_out = fname_fem_out.substr(0,fname_mat.length()-4);
	fname_fem_out += "_fem_reordered.bin";
	std::ofstream os (fname_fem_out.c_str(), std::ios::binary);
	CS_write(os, n, mat.getRowPtrArray(), mat.getColIdxArray(), mat.getEntryArray());
#endif

	if (verbose) {
		std::cout << "matrix vector multiplication with Toms amuxCRS ... " << std::flush;
	}
	run_timer.start();
	cpu_timer.start();
	for (size_t k(0); k<n_mults; k++) {
		mat.amux (1.0, x, y);
	}
	cpu_timer.stop();
	run_timer.stop();

	if (verbose) {
		std::cout << "done [" << cpu_timer.elapsed() << " sec cpu time], ["
				<< run_timer.elapsed() << " sec run time]" << std::endl;
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
