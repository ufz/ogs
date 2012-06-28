/**
 * MatVecMultNDPermOpenMP.cpp
 *
 * Created on 2012-01-20 by Thomas Fischer
 */

#include <cstdlib>

// BaseLib
#include "RunTime.h"
#include "CPUTime.h"

// MathLib
#include "sparse.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReorderedOpenMP.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"

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
		BaseLib::RunTime timer;
		timer.start();
		CS_read(in, n, iA, jA, A);
		timer.stop();
		if (verbose) {
			std::cout << "ok, [wclock: " << timer.elapsed() << " s]" << std::endl;
		}
	} else {
		std::cout << "error reading matrix from " << fname_mat << std::endl;
	}
	unsigned nnz(iA[n]);
	if (verbose) {
		std::cout << "Parameters read: n=" << n << ", nnz=" << nnz << std::endl;
	}

	MathLib::CRSMatrixReorderedOpenMP mat(n, iA, jA, A);

	double *x(new double[n]);
	double *y(new double[n]);

	for (unsigned k(0); k<n; ++k)
		x[k] = 1.0;

	// create time measurement objects
	BaseLib::RunTime run_timer;
	BaseLib::CPUTime cpu_timer;

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
	if (verbose) {
		std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
	} else {
		if (argc == 4) {
			std::ofstream result_os(argv[3], std::ios::app);
			if (result_os) {
				result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << " calc nested dissection perm" << std::endl;
			}
			result_os.close();
		} else {
			std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << " calc nested dissection perm" << std::endl;
		}
	}


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
	else {
		if (argc == 4) {
			std::ofstream result_os(argv[3], std::ios::app);
			if (result_os) {
				result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << " applying nested dissection perm" << std::endl;
			}
			result_os.close();
		} else {
			std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
		}
	}

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
		std::cout << "done [" << cpu_timer.elapsed() << " sec cpu time], [wclock: "
				<< run_timer.elapsed() << " sec]" << std::endl;
	} else {
		if (argc == 4) {
			std::ofstream result_os (argv[3], std::ios::app);
			if (result_os) {
				result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << " " << n_mults << " MatVecMults, matrix " << fname_mat << std::endl;
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
