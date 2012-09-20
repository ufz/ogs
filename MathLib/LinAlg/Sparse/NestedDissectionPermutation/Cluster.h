/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file Cluster.h
 *
 * Created on 2012-01-02 by Thomas Fischer
 */

#ifndef CLUSTER_H_
#define CLUSTER_H_

#include "ClusterBase.h"

namespace MathLib {

/** \brief class for storing clusters of degrees of freedom without
    geometric information. This class stores clusters (not separators).
*/

class Cluster: public ClusterBase
{
public:
	/**
	 * Constructor creates the root of the cluster tree
	 * @param n
	 * @param jA
	 * @param iA
	 * @return
	 */
	Cluster(unsigned n, unsigned* iA, unsigned* jA);

	virtual void subdivide(unsigned bmin);

	/** Method returns the status of this ClusterBase object. In this case
	 * instances of this class are "normal" Clusters.
	 * @return false
	 */
	virtual bool isSeparator() const
	{
		return false;
	}

	/** Destructor. */
	virtual ~Cluster() {}

	/**
	 * Method creates recursively the cluster tree, i.e. changes the permutation
	 * op_perm and po_perm and create child cluster trees. For this task only the
	 * adjacency matrix is used.
	 * @param op_perm permutation: original_idx = op_perm[permutated_idx]
	 * @param po_perm reverse permutation: permutated_idx = po_perm[original_idx]
	 * @param bmin threshold value for stopping further refinement
	 * @return a cluster tree
	 */
	virtual void createClusterTree(unsigned* op_perm, unsigned* po_perm,
			unsigned bmin = 50);

protected:
	/** \brief Constructor
	 \param father parent node in cluster tree
	 \param beg beginning index of the cluster
	 \param end beginning index of the next cluster
	 \param op_perm permutation
	 \param po_perm permutation
	 \param global_mat reference to adjacency matrix of the matrix graph in
	 crs format
	 \param local_mat pointer to the local adjacency matrix of the matrix
	 graph in crs format
	 */
	Cluster(ClusterBase* father, unsigned beg, unsigned end, unsigned* op_perm,
			unsigned* po_perm, AdjMat* global_mat, AdjMat* local_mat);

private:
	/** update perm */
	void updatePerm(unsigned* reordering, unsigned &isep0, unsigned &isep1, unsigned* l_op_perm, unsigned* l_po_perm);
};

} // end namespace MathLib

#endif /* CLUSTER_H_ */
