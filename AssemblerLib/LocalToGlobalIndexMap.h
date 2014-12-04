/**
 * \file LocalToGlobalIndexMap.h
 * \author Norihiro Watanabe
 * \author Wenqing Wang
 * \date   2013-04-16, 2014-11-14
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
#define ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_

#include <vector>

#ifdef USE_PETSC
#include <petscsys.h>
#endif

#include "AssemblerLib/MeshComponentMap.h"
#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MeshLib/MeshSubsets.h"

namespace AssemblerLib
{

/// Row and column indices in global linear algebra objects for each mesh item.
///
/// The row and column indices in global matrix and rhs vector for example,
/// required for addition of local contributions from every mesh item (like node
/// or cell) to global objects.
///
/// The number of rows should be equal to the number of mesh items and the
/// number of columns should be equal to the number of the components on that
/// mesh item.
class LocalToGlobalIndexMap
{
public:
#ifdef USE_PETSC
    typedef MathLib::RowColumnIndices<PetscInt> RowColumnIndices;
#else
    typedef MathLib::RowColumnIndices<std::size_t> RowColumnIndices;
#endif    
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubsets*> const& mesh_subsets,
        AssemblerLib::ComponentOrder const order =
            AssemblerLib::ComponentOrder::BY_COMPONENT,
            const bool is_linear_element = true);

    /// Returns total number of degrees of freedom.
    std::size_t dofSize() const;

    /// Returns total number of global degrees of freedom for DDC.
    std::size_t dofSizeGlobal() const
    {
		return _mesh_component_map.getGlobalDOF();
    }

    std::size_t size() const;

    RowColumnIndices operator[](std::size_t const mesh_item_id) const;

    LineIndex rowIndices(std::size_t const mesh_item_id) const;
    LineIndex columnIndices(std::size_t const mesh_item_id) const;
    
#ifdef USE_PETSC
    const std::vector<bool> &getNodeGhostFlags(std::size_t const mesh_item_id) const
    {
        return _element_ghost_node_flags[mesh_item_id];
    }
#endif    

private:
    std::vector<MeshLib::MeshSubsets*> const& _mesh_subsets;
    AssemblerLib::MeshComponentMap _mesh_component_map;

    /// Vector contains for each element a vector of global row/or entry indices
    /// in the global stiffness matrix or vector
    std::vector<LineIndex> _rows;

    /// Vector alias to that contains for each element a vector of global column indices
    /// in the global stiffness matrix
    const std::vector<LineIndex> &_columns = _rows;

#ifdef USE_PETSC
    /// Element nodal flag for whether the node is a ghost node.
    /// Its first entry of each sub-vector is the indicator for whether the element is a ghost element
    std::vector<std::vector<bool>> _element_ghost_node_flags;
#endif
    
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
