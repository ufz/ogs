/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Separator.h
 *
 * Created on 2012-01-02 by Thomas Fischer
 */

#ifndef SEPARATOR_H_
#define SEPARATOR_H_

#include "LinAlg/Sparse/NestedDissectionPermutation/ClusterBase.h"

namespace MathLib {

class Cluster;
class AdjMat;

/** \brief class for storing clusters of degrees of freedom without
 geometric information.
*/
class Separator: public ClusterBase
{
public:
	/** brief Constructor builds a initial object for clustering
	 \param father pointer to the father node in cluster tree
	 \param beg index in op_perm and po_perm which describes the begin of the index set of the Separator
	 \param end index in op_perm and po_perm which describes the begin of the index set of the next
	 ClusterBase
	 \param op_perm permutation
	 \param po_perm inverse permutation
	 \param global_mat reference to adjacency matrix of the matrix graph in crs format
	 \param local_mat reference to the local adjacency matrix of the matrix graph in crs format
	 */
	Separator(ClusterBase* father, unsigned beg, unsigned end,
			unsigned* op_perm, unsigned* po_perm, AdjMat* global_mat,
			AdjMat* local_mat);

	/** Destructor. */
	virtual ~Separator();

	/** Method returns the status of this ClusterAlg object. Instances
	 of this class are Separators.
	 \returns true
	 */
	virtual bool isSeparator() const
	{
		return true;
	}

private:
	/** update perm */
	unsigned updatePerm(unsigned *reordering, unsigned* l_op_perm, unsigned* l_po_perm);
};

} // end namespace MathLib

#endif /* SEPARATOR_H_ */
