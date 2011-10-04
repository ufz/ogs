#include "sparse.h"

namespace MathLib {

bool generateDiagPrecond (unsigned n, unsigned* iA, unsigned* jA, double* A,
	double* diag)
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

} // end namespace MathLib
