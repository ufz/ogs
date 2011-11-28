#include "sparse.h"
#include <limits>
#include <cmath>

namespace MathLib {

bool generateDiagPrecond (unsigned n, unsigned const*const iA, unsigned const*const jA,
				double const*const A, double* diag)
{
	unsigned idx; // first idx of next row
	unsigned c; // column
	unsigned j;
	bool has_no_diag;

	for (unsigned r(0); r<n; ++r) {
		idx=iA[r+1];
		has_no_diag=true;
		for (j=iA[r]; j<idx && has_no_diag; ++j) {
			c=jA[j];
			if (c==r) {
				has_no_diag=false;
				diag[r] = 1.0/A[j];
			}
		}
		if (j==idx && has_no_diag) {
			std::cout << "row " << r << " has no diagonal element " << std::endl;
			return false;
		}
	}
	return true;
}

bool generateDiagPrecondRowSum(unsigned n, unsigned const*const iA, unsigned const*const jA,
				double const*const A, double* diag)
{
	unsigned idx; // first idx of next row
	unsigned c; // column
	unsigned j;

	for (unsigned r(0); r<n; ++r) {
		diag[r] = 0.0;
		idx=iA[r+1];
		for (j=iA[r]; j<idx; ++j) {
			diag[r] += fabs(A[j]);
		}
		if (fabs(diag[r]) < std::numeric_limits<double>::epsilon()) {
			std::cout << "row " << r << " has only very small entries" << std::endl;
			return false;
		}
		diag[r] = 1.0/diag[r];
	}
	return true;
}

bool generateDiagPrecondRowMax(unsigned n, unsigned const*const iA, unsigned const*const jA,
				double const*const A, double* diag)
{
	unsigned idx; // first idx of next row
	unsigned c; // column
	unsigned j;

	for (unsigned r(0); r<n; ++r) {
		idx=iA[r+1];
		diag[r] = A[idx];
		for (j=iA[r]; j<idx; ++j) {
			if (A[j] > diag[r])
				diag[r] = A[j];
		}
		if (fabs(diag[r]) < std::numeric_limits<double>::epsilon()) {
			std::cout << "the maximum entry of row " << r << " has only very small value" << std::endl;
			return false;
		}
		diag[r] = 1.0/diag[r];
	}
	return true;
}

} // end namespace MathLib
