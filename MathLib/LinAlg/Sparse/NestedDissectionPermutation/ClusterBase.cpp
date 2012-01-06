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
	_g_po_perm(NULL), _g_adj_mat(NULL), _l_adj_mat(NULL)
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

	// make a copy of the local row_ptr array
	unsigned const* l_row_ptr(_l_adj_mat->getRowPtrArray());
	unsigned *g_row_ptr(new unsigned[n + 1]);
	for (unsigned k = 0; k <= n; ++k)
		g_row_ptr[k] = l_row_ptr[k];
	// make a copy of the local col_idx array
	unsigned const* l_col_idx(_l_adj_mat->getColIdxArray());
	const unsigned g_nnz(g_row_ptr[n]);
	unsigned *g_col_idx(new unsigned[g_nnz]);
	for (unsigned k = 0; k < g_nnz; ++k)
		g_col_idx[k] = l_col_idx[k];
	// generate global matrix from local matrix
	// (only in the root of cluster tree)
	_g_adj_mat = new AdjMat(n, g_row_ptr, g_col_idx);


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
