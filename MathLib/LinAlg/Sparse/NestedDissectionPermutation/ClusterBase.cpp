/*
 * ClusterBase.cpp
 *
 *  Created on: 02.01.2012
 *      Author: TF
 */

//#include "blas.h"
#include "AdjMat.h"
#include "LinAlg/Sparse/NestedDissectionPermutation/ClusterBase.h"

namespace MathLib {

ClusterBase::ClusterBase(unsigned n, unsigned const* const iA,
		unsigned const*const jA) :
	_beg(0), _end(n), _n_sons(0), _sons(NULL), _parent(NULL), _g_op_perm(NULL),
	_g_po_perm(NULL), _l_adj_mat(NULL)
{
	const unsigned nnz = iA[n];

	// create adjacency matrix
	unsigned *row_ptr = new unsigned[n + 1];
	for (unsigned k = 0; k <= n; ++k)
		row_ptr[k] = iA[k];
	unsigned *col_idx = new unsigned[nnz];
	for (unsigned k = 0; k < nnz; ++k)
		col_idx[k] = jA[k];

	_l_adj_mat = new AdjMat(n, row_ptr, col_idx);
	_l_adj_mat->makeSymmetric();
}

ClusterBase::ClusterBase(ClusterBase *father, unsigned beg, unsigned end,
		unsigned* op_perm, unsigned* po_perm, AdjMat* global_mat, AdjMat* local_mat) :
	_beg(beg), _end(end), _n_sons(0), _sons(NULL), _parent(father),
	_g_op_perm(op_perm), _g_po_perm(po_perm), _g_adj_mat(global_mat),
	_l_adj_mat(local_mat)
{
}

ClusterBase::~ClusterBase()
{
	if (_parent == NULL)
		delete _g_adj_mat;
	delete _l_adj_mat;
}

} // end namespace MathLib
