#ifndef CRSMATRIXDIAGPRECOND_H
#define CRSMATRIXDIAGPRECOND_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CRSMatrix.h"
#include "sparse.h"
#include "../Preconditioner/generateDiagPrecond.h"

namespace MathLib {

/**
 * Class CRSMatrixDiagPrecond represents a matrix in compressed row storage
 * format associated with a diagonal preconditioner.
 *
 * The user can either read the matrix from a file or give corresponding arrays
 * to an alternative constructor. In both cases the user have to calculate the
 * preconditioner explicit via calcPrecond() method!
 */
class CRSMatrixDiagPrecond : public CRSMatrix<double, unsigned>
{
public:
	/**
	 * Constructor takes a file name. The file is read in binary format
	 * by the constructor of the base class (template) CRSMatrix.
	 * Further details you can see in function CS_read().
	 *
	 * The user have to calculate the preconditioner explicit via calcPrecond() method!
	 *
	 * @param fname the name of the file that contains the matrix in
	 * binary compressed row storage format
	 */
	CRSMatrixDiagPrecond(std::string const &fname) :
		CRSMatrix<double, unsigned> (fname), _inv_diag(NULL)
	{}

	/**
	 * Constructs a matrix object from given data.
	 *
	 * The user have to calculate the preconditioner explicit via calcPrecond() method!
	 * @param n number of rows / columns of the matrix
	 * @param iA row pointer of matrix in compressed row storage format
	 * @param jA column index of matrix in compressed row storage format
	 * @param A data entries of matrix in compressed row storage format
	 */
	CRSMatrixDiagPrecond(unsigned n, unsigned *iA, unsigned *jA, double* A) :
		CRSMatrix<double, unsigned> (n, iA, jA, A), _inv_diag(NULL)
	{}

	void calcPrecond()
	{
		if (_inv_diag != NULL)
			delete [] _inv_diag;
		_inv_diag = new double[_n_rows];

		if (!generateDiagPrecond(_n_rows, _row_ptr, _col_idx, _data, _inv_diag)) {
			std::cout << "Could not create diagonal preconditioner" << std::endl;
		}
//		if (!generateDiagPrecondRowSum(_n_rows, _row_ptr, _data, _inv_diag)) {
//			std::cout << "Could not create diagonal preconditioner" << std::endl;
//		}
//		if (!generateDiagPrecondRowMax(_n_rows, _row_ptr, _data, _inv_diag)) {
//			std::cout << "Could not create diagonal preconditioner" << std::endl;
//		}

	}

	void precondApply(double* x) const
	{
		for (unsigned k=0; k<_n_rows; ++k) {
			x[k] = _inv_diag[k]*x[k];
		}
	}

	~CRSMatrixDiagPrecond()
	{
		delete [] _inv_diag;
	}
private:
	double *_inv_diag;
};

}

#endif

