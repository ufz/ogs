#ifndef SPARSEMATRIXBASE_H
#define SPARSEMATRIXBASE_H

template<class T> class SparseMatrixBase
{
public:
	SparseMatrixBase(unsigned n1, unsigned n2) : _n_rows(n1), _n_cols(n2) { }
	SparseMatrixBase() : _n_rows(0), _n_cols(0) { }
	virtual void amux(T d, T const * const x, T *y) const = 0;         // y +=d*Ax
	virtual ~SparseMatrixBase() { }
	unsigned getNRows () const { return _n_rows; }
	unsigned getNCols () const { return _n_cols; }

protected:
	unsigned _n_rows; 
	unsigned _n_cols; 
};


#endif

