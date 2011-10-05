#ifndef CRSMATRIXDIAGPRECOND_H
#define CRSMATRIXDIAGPRECOND_H

#include <omp.h>

#include "CRSMatrix.h"
#include "sparse.h"
#include "../Preconditioner/generateDiagPrecond.h"

namespace MathLib {

class CRSMatrixDiagPrecond : public CRSMatrix<double>
{
public:
	CRSMatrixDiagPrecond (std::string const &fname)
		: CRSMatrix<double>(fname), _inv_diag(new double[_n_rows])
	{
		if (!generateDiagPrecond (_n_rows, _row_ptr, _col_idx, _data, _inv_diag)) {
			std::cout << "Could not create diagonal preconditioner" << std::endl;
		}
	}

	void precondApply(double* x) const {
		{
			unsigned k;
			#pragma omp parallel for
			for (k=0; k<_n_rows; ++k) {
				x[k] = _inv_diag[k]*x[k];
			}
		}
	}

	~CRSMatrixDiagPrecond() {
		delete [] _inv_diag;
	}
private:
	double *_inv_diag;
};

}

#endif

