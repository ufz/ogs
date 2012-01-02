/*
 * Separator.cpp
 *
 *  Created on: 02.01.2012
 *      Author: TF
 */

// BaseLib
#include "swap.h"

#include "LinAlg/Sparse/NestedDissectionPermutation/Separator.h"

namespace MathLib {

extern "C" void METIS_PartGraphRecursive(unsigned*, unsigned*, unsigned*,
					 unsigned*, unsigned*, unsigned*,
					 unsigned*, unsigned*, unsigned*,
					 unsigned*, unsigned*);

Separator::Separator(ClusterBase* father, unsigned beg, unsigned end,
                     unsigned* op_perm, unsigned* po_perm,
                     AdjMat* global_mat, AdjMat* local_mat)
  : ClusterBase (father, beg, end, op_perm, po_perm, global_mat, local_mat)
{}

Separator::~Separator()
{}

unsigned Separator::updatePerm(unsigned* reordering, unsigned* l_op_perm, unsigned* l_po_perm)
{
  unsigned beg = 0, end = _end-_beg;
  while (beg < end) {
    if (reordering[beg] == 1) {
      --end;
      while (beg < end && reordering[end] == 1) --end;
      // local permutation
      BaseLib::swap(l_op_perm [beg], l_op_perm[end]);
      BaseLib::swap(l_po_perm[l_op_perm [beg]], l_po_perm[l_op_perm[end]]);
    }
    ++beg;
  }
  return ((beg>end) ? beg-1 : end);
}

} // end namespace MathLib
