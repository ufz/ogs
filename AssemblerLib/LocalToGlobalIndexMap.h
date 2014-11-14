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
    typedef MathLib::RowColumnIndices<std::size_t> RowColumnIndices;
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

    std::size_t size() const;

    RowColumnIndices operator[](std::size_t const mesh_item_id) const;

    LineIndex rowIndices(std::size_t const mesh_item_id) const;
    LineIndex columnIndices(std::size_t const mesh_item_id) const;

private:
    std::vector<MeshLib::MeshSubsets*> const& _mesh_subsets;
    AssemblerLib::MeshComponentMap _mesh_component_map;

    /// Vector contains for each element a vector of global row/or entry indices
    /// in the global stiffness matrix or vector
    std::vector<LineIndex> _rows;

    /// Vector contains for each element a vector of global column indices
    /// in the global stiffness matrix
    std::vector<LineIndex> _columns_real;

    /// Vector alias to that contains for each element a vector of global column indices
    /// in the global stiffness matrix
    std::vector<LineIndex> &_columns = _rows;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
