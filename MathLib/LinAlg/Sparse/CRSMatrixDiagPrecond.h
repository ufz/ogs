#ifndef CRSMATRIXDIAGPRECOND_H
#define CRSMATRIXDIAGPRECOND_H

#include <omp.h>

#include "CRSMatrix.h"
#include "sparse.h"

class CRSMatrixDiagPrecond : public CRSMatrix<double>
{
public:
	CRSMatrixDiagPrecond (std::string const &fname, unsigned num_of_threads=1) 
		: CRSMatrix<double>(fname, num_of_threads), _inv_diag(NULL) 
	{
		_inv_diag = new double[_n_rows];
		if (!generateDiagPrecond (_n_rows, _row_ptr, _col_idx, _data, _inv_diag)) {
			std::cout << "Could not create diagonal preconditioner" << std::endl;
		}
	}

	void precondApply(double* x) const {
		omp_set_num_threads( _num_of_threads );
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

#endif

