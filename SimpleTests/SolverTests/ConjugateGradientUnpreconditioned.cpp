#include <fstream>
#include <iostream>
#include <omp.h>
#include "LinAlg/Solvers/CG.h"
#include "LinAlg/Sparse/CRSMatrix.h"
#include "sparse.h"
#include "vector_io.h"
//#include "timeMeasurement.h"

int main(int argc, char *argv[])
{
	// *** reading matrix in crs format from file
	std::string fname("/work/fischeth/data/testmat.bin");
//	std::ifstream in(fname.c_str(), std::ios::binary);
	MathLib::CRSMatrix<double, unsigned> *mat (new MathLib::CRSMatrix<double, unsigned>(fname));
/*	double *A(NULL);
	unsigned *iA(NULL), *jA(NULL), n;
	if (in) {
		std::cout << "reading matrix from " << fname << " ... " << std::flush;
		CS_read(in, n, iA, jA, A);
		in.close();
		std::cout << " ok" << std::endl;
	}
	unsigned nnz(iA[n]);
*/
	unsigned n (mat->getNRows());
	std::cout << "Parameters read: n=" << n << std::endl;

	double *x(new double[n]);
	double *b(new double[n]);

	// *** init start vector x
	for (size_t k(0); k<n; k++) {
		x[k] = 1.0;
	}
	// *** read rhs
	fname = "/work/fischeth/data/rhs.dat";
	std::ifstream in(fname.c_str());
	if (in) {
		read (in, n, b);
		in.close();
	} else {
		std::cout << "problem reading rhs - initializing b with 1.0" << std::endl;
		for (size_t k(0); k<n; k++) {
			b[k] = 1.0;
		}
	}

	std::cout << "solving system with CG method ... " << std::flush;
	time_t start_time, end_time;
	time(&start_time);
//	double cg_time (cputime(0.0));
	double eps (1.0e-6);
	unsigned steps (4000);
	CG (mat, b, x, eps, steps, 1);
//	cg_time = cputime(cg_time);
	time(&end_time);
	std::cout << " in " << steps << " iterations (residuum is " << eps << ") took " << /*cg_time <<*/ " sec time and " << (end_time-start_time) << " sec" << std::endl;

	delete mat;
	delete [] x;
	delete [] b;

	return 0;
}

