#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE> class SparseMatrixBase
{
public:
	SparseMatrixBase(IDX_TYPE n1, IDX_TYPE n2) :
			_n_rows(n1), _n_cols(n2)
	{}
	SparseMatrixBase() :
			_n_rows(static_cast<IDX_TYPE>(0)), _n_cols(static_cast<IDX_TYPE>(0))
	{}
	/**
	 * y = d * A * x
	 * @param d scalar factor
	 * @param x vector to multiply with
	 * @param y result vector
	 */
	virtual void amux(FP_TYPE d, FP_TYPE const * const __restrict__ x, FP_TYPE * __restrict__ y) const = 0;
	virtual ~SparseMatrixBase() {}
	/**
	 * get the number of rows
	 * @return the number of rows
	 */
	IDX_TYPE getNRows () const { return _n_rows; }
	/**
	 * get the number of columns
	 * @return the number of columns
	 */
	IDX_TYPE getNCols () const { return _n_cols; }

protected:
	/**
	 * the number of rows
	 */
	IDX_TYPE _n_rows;
	/**
	 * the number of columns
	 */
	IDX_TYPE _n_cols;
};

} // end namespace MathLib

#endif
