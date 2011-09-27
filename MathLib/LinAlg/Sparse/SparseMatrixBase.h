#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

#include "../MatrixBase.h"

namespace MathLib {

template<class T> class SparseMatrixBase : public MatrixBase
{
public:
	SparseMatrixBase(unsigned n1, unsigned n2) : MatrixBase (n1,n2) {}
	SparseMatrixBase() : MatrixBase () {}
	virtual void amux(T d, T const * const x, T *y) const = 0;         // y +=d*Ax
	virtual ~SparseMatrixBase() { }
};

} // end namespace MathLib

#endif
