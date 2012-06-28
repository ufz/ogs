/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file ClusterBase.h
 *
 * Created on 2012-01-02 by Thomas Fischer
 */

#ifndef CLUSTERBASE_H_
#define CLUSTERBASE_H_

namespace MathLib {

class AdjMat;

/** \brief Base class for storing clusters of degrees of freedom without
 geometric information.
 */
class ClusterBase
{
public:
	/**
	 * Constructor creates the root of the cluster tree
	 * @param n number of rows/columns
	 * @param iA row pointer array
	 * @param jA column index array
	 * @return
	 */
	ClusterBase(unsigned n, unsigned const*const iA, unsigned const*const jA);
	/*!
	 \brief Constructor
	 \param father parent node in cluster tree
	 \param beg beginning index of the cluster
	 \param end beginning index of the next cluster
	 \param op_perm global permutation array (original_idx = op_perm[permuted_idx])
	 \param po_perm global permutation array (permuted_idx = po_perm[original_idx])
	 \param global_mat reference to the global adjacency matrix of the matrix graph in crs format
	 \param local_mat pointer to the local adjacency matrix of the matrix graph in crs format
	 */
	ClusterBase(ClusterBase* father, unsigned beg, unsigned end,
			unsigned* op_perm, unsigned* po_perm, AdjMat* global_mat, AdjMat* local_mat);

	/** \brief Destructor.
	 * Destructor frees all form the objects allocated memory.
	 * */
	virtual ~ClusterBase();

	virtual bool isSeparator() const = 0;

#ifndef NDEBUG
	AdjMat const* getGlobalAdjMat() const { return _g_adj_mat; }
#endif

protected:
	/** \brief Method returns the pointer to the parent cluster.
	 \returns parent cluster */
	ClusterBase* getParent() const
	{
		return _parent;
	}
	/**
	 * beginning index in the global permutation arrays
	 */
	const unsigned _beg;

	/**
	 * beginning index of next next cluster in the global permutation arrays
	 */
	const unsigned _end;

	/**
	 * number of sons, _n_sons==0 iff this is a leaf
	 */
	unsigned _n_sons;

	/**
	 * the array of sons of this cluster in the cluster tree
	 */
	ClusterBase** _sons;

	/**
	 * pointer to parent
	 */
	ClusterBase *_parent;
	/**
	 * pointer global permutation array (original_idx = op_perm[permuted_idx])
	 */
	unsigned* _g_op_perm;
	/**
	 * global permutation: permutation <- po_perm <- original
	 */
	unsigned* _g_po_perm;
	/**
	 * pointer to global adjacency matrix
	 * The attribute _g_adj_mat stores the set of edges of the matrix graph $G = (V,E)$.
	 * (see class AdjMat)
	 */
	AdjMat* _g_adj_mat;
	/**
	 * local adjacency matrix
	 */
	AdjMat* _l_adj_mat;
};

} // end namespace MathLib

#endif /* CLUSTERBASE_H_ */
