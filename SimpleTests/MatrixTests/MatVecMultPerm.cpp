/*
 * MatVecMultPerm.cpp
 *
 *  Created on: Jan 3, 2012
 *      Author: TF
 */

#include <cstdlib>

// BaseLib
#include "RunTime.h"
#include "CPUTime.h"
#include "tclap/CmdLine.h"

// MathLib
#include "sparse.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/CRSMatrixReordered.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"
#include "LinAlg/Sparse/CRSMatrix.h"

int main(int argc, char *argv[])
{
	TCLAP::CmdLine cmd("Simple matrix vector multiplication test", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m matrix".
	TCLAP::ValueArg<std::string> matrix_arg("m","matrix","input matrix file",true,"","string");

	// Add the argument mesh_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( matrix_arg );

//	TCLAP::ValueArg<unsigned> n_cores_arg("n", "number-cores", "number of cores to use", true, "1", "number");
//	cmd.add( n_cores_arg );

	TCLAP::ValueArg<unsigned> n_mults_arg("n", "number-of-multiplications", "number of multiplications to perform", true, 10, "unsigned");
	cmd.add( n_mults_arg );

	TCLAP::ValueArg<std::string> output_arg("o", "output", "output file", false, "", "string");
	cmd.add( output_arg );

	TCLAP::ValueArg<unsigned> verbosity_arg("v", "verbose", "level of verbosity [0 very low information, 1 much information]", false, 0, "string");
	cmd.add( verbosity_arg );

	cmd.parse( argc, argv );

	// read the number of multiplication to execute
	unsigned n_mults (n_mults_arg.getValue());
	std::string fname_mat (matrix_arg.getValue());
	bool verbose (verbosity_arg.getValue());

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

//	MathLib::CRSMatrix<double, unsigned> mat(n, iA, jA, A);
	MathLib::CRSMatrixReordered mat(n, iA, jA, A);

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
		if (! output_arg.getValue().empty()) {
			std::ofstream result_os(output_arg.getValue().c_str(), std::ios::app);
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
		if (! ((output_arg.getValue()).empty())) {
			std::ofstream result_os((output_arg.getValue()).c_str(), std::ios::app);
			if (result_os) {
				result_os << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << " applying nested dissection perm" << std::endl;
			}
			result_os.close();
		} else {
			std::cout << cpu_timer.elapsed() << "\t" << run_timer.elapsed() << std::endl;
		}
	}

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
//	std::string fname_fem_out (fname_mat);
//	fname_fem_out = fname_fem_out.substr(0,fname_mat.length()-4);
//	fname_fem_out += "_fem_reordered.bin";
//	std::ofstream os (fname_fem_out.c_str(), std::ios::binary);
//	CS_write(os, n, mat.getRowPtrArray(), mat.getColIdxArray(), mat.getEntryArray());
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
		std::cout << "done [" << cpu_timer.elapsed() << " sec cpu time], [wclock: "
				<< run_timer.elapsed() << " sec]" << std::endl;
	} else {
		if (! output_arg.getValue().empty()) {
			std::ofstream result_os (output_arg.getValue().c_str(), std::ios::app);
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
