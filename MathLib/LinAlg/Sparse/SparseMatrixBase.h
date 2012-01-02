#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

#include "../MatrixBase.h"

namespace MathLib {

template<typename FP_TYPE, typename IDX_TYPE> class SparseMatrixBase : public MatrixBase
{
public:
	SparseMatrixBase(IDX_TYPE n1, IDX_TYPE n2) : MatrixBase (n1,n2) {}
	SparseMatrixBase() : MatrixBase () {}
	/**
	 * y = d * A * x
	 * @param d scalar factor
	 * @param x vector to multiply with
	 * @param y result vector
	 */
	virtual void amux(FP_TYPE d, FP_TYPE const * const __restrict__ x, FP_TYPE * __restrict__ y) const = 0;
	virtual ~SparseMatrixBase() {};
};

} // end namespace MathLib

#endif
