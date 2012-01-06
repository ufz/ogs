/*
 * Cluster.cpp
 *
 *  Created on: 02.01.2012
 *      Author: TF
 */

#include "metis.h"

// BaseLib
#include "swap.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/Cluster.h"
//#include "blas.h"
#include "Cluster.h"
#include "Separator.h"
#include "AdjMat.h"

// METIS function
//extern "C" void METIS_ComputeVertexSeparator(unsigned*, unsigned*, unsigned*,
//					   unsigned*, unsigned*, unsigned*, unsigned*);

//extern "C" void METIS_NodeND(unsigned*, unsigned*, unsigned*,
//				   unsigned*, unsigned*, unsigned*, unsigned*);

namespace MathLib {

Cluster::Cluster (unsigned n, unsigned* iA, unsigned* jA)
  : ClusterBase (n, iA, jA)
{}


Cluster::Cluster(ClusterBase* father, unsigned beg, unsigned end,
                         unsigned* op_perm, unsigned* po_perm,
                         AdjMat* global_mat, AdjMat* local_mat)
  : ClusterBase(father, beg, end, op_perm, po_perm, global_mat, local_mat)
{}

void Cluster::subdivide(unsigned bmin)
{
	const unsigned size(_end - _beg);
	if (size > bmin) {

		unsigned n_rows(_l_adj_mat->getNRows());
		unsigned *xadj(const_cast<unsigned*>(_l_adj_mat->getRowPtrArray()));
		unsigned *adjncy(const_cast<unsigned*>(_l_adj_mat->getColIdxArray()));
//		unsigned nparts = 2;
		unsigned options[METIS_NOPTIONS]; // for METIS
		METIS_SetDefaultOptions(options);
//		options[METIS OPTION PTYPE] = METIS PTYPE RB;
//		options[METIS OPTION OBJTYPE] = METIS OBJTYPE CUT;
//		options[METIS OPTION CTYPE] = METIS CTYPE SHEM;
//		options[] = ;
//		options[] = ;
//		options[] = ;

//		unsigned sepsize(0); // for METIS
		unsigned *vwgt(new unsigned[n_rows + 1]);
//		const unsigned nnz(xadj[n_rows]);
//		unsigned *adjwgt(new unsigned[nnz]);
		for (unsigned k(0); k < n_rows + 1; k++)
			vwgt[k] = 1;
//		for (unsigned k(0); k < nnz; k++)
//			adjwgt[k] = 1;
//		unsigned *part(new unsigned[n_rows + 1]);

		// subdivide the index set into three parts employing METIS
//		METIS_ComputeVertexSeparator(&n_rows, xadj, adjncy, vwgt, &options,
//				&sepsize, part);
		METIS_NodeND(&n_rows, xadj, adjncy, vwgt, options, _g_op_perm, _g_po_perm);

//		// create and init local permutations
//		unsigned *l_op_perm(new unsigned[size]);
//		unsigned *l_po_perm(new unsigned[size]);
//		for (unsigned i = 0; i < size; ++i)
//			l_op_perm[i] = l_po_perm[i] = i;
//
//		unsigned isep1, isep2;
//		updatePerm(part, isep1, isep2, l_op_perm, l_po_perm);
//		delete[] part;
//
//		// update global permutation
//		unsigned *t_op_perm = new unsigned[size];
//		for (unsigned k = 0; k < size; ++k)
//			t_op_perm[k] = _g_op_perm[_beg + l_op_perm[k]];
//
//		for (unsigned k = _beg; k < _end; ++k) {
//			_g_op_perm[k] = t_op_perm[k - _beg];
//			_g_po_perm[_g_op_perm[k]] = k;
//		}
//		delete[] t_op_perm;
//
//		// next recursion step
//		if ((isep1 >= bmin) && (isep2 - isep1 >= bmin)) {
//			// construct adj matrices for [0, isep1), [isep1,isep2), [isep2, _end)
//			AdjMat *l_adj0(_l_adj_mat->getMat(0, isep1, l_op_perm, l_po_perm));
//			AdjMat *l_adj1(_l_adj_mat->getMat(isep1, isep2, l_op_perm, l_po_perm));
//			AdjMat *l_adj2(_l_adj_mat->getMat(isep2, size, l_op_perm, l_po_perm));
//
//			delete[] l_op_perm;
//			delete[] l_po_perm;
//			delete _l_adj_mat;
//			_l_adj_mat = NULL;
//
//			_n_sons = 3;
//			_sons = new ClusterBase*[_n_sons];
//
//			isep1 += _beg;
//			isep2 += _beg;
//
//			// constructing child nodes for index cluster tree
//			_sons[0] = new Cluster(this, _beg, isep1, _g_op_perm, _g_po_perm, _g_adj_mat, l_adj0);
//			_sons[1] = new Cluster(this, isep1, isep2, _g_op_perm, _g_po_perm, _g_adj_mat, l_adj1);
//			_sons[2] = new Separator(this, isep2, _end, _g_op_perm,	_g_po_perm, _g_adj_mat, l_adj2);
//
//			dynamic_cast<Cluster*>(_sons[0])->subdivide(bmin);
//			dynamic_cast<Cluster*>(_sons[1])->subdivide(bmin);
//
//		} else {
//			delete _l_adj_mat;
//			_l_adj_mat = NULL;
//		} // end if next recursion step
	} // end if ( connected && size () > bmin )
}


void Cluster::updatePerm(unsigned* reordering, unsigned &isep0,
		unsigned &isep1, unsigned* l_op_perm, unsigned* l_po_perm)
{
	unsigned beg = 0, end = _end - _beg;
	while (beg < end) {
		if (reordering[beg] >= 1) {
			--end;
			while (beg < end && reordering[end] >= 1)
				--end;
			// local permutation
			BASELIB::swap(l_op_perm[beg], l_op_perm[end]);
			BASELIB::swap(l_po_perm[l_op_perm[beg]], l_po_perm[l_op_perm[end]]);
			BASELIB::swap(reordering[beg], reordering[end]);
		}
		++beg;
	}
	if (beg > end)
		isep0 = beg - 1;
	else
		isep0 = end;

	beg = isep0, end = _end - _beg;
	while (beg < end) {
		if (reordering[beg] == 2) {
			--end;
			while (beg < end && reordering[end] == 2)
				--end;
			// local permutation
			BASELIB::swap(l_op_perm[beg], l_op_perm[end]);
			BASELIB::swap(l_po_perm[l_op_perm[beg]], l_po_perm[l_op_perm[end]]);
			BASELIB::swap(reordering[beg], reordering[end]);
		}
		++beg;
	}
	if (beg > end)
		isep1 = beg - 1;
	else
		isep1 = end;
}


void Cluster::createClusterTree(unsigned* op_perm, unsigned* po_perm,
		unsigned bmin)
{
	_g_op_perm = op_perm;
	_g_po_perm = po_perm;

	// *** 1 create local problem
	unsigned n = _g_adj_mat->getNRows();
	unsigned *l_op_perm = new unsigned[n];
	unsigned *l_po_perm = new unsigned[n];
	for (unsigned k = 0; k < n; ++k)
		l_op_perm[k] = l_po_perm[k] = k;
	_l_adj_mat = _l_adj_mat->getMat(0, n, l_op_perm, l_po_perm);

	// *** 2 create cluster tree
	subdivide(bmin);
}

} // end namespace MathLib
