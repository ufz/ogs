/**
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
#define ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_

#include <vector>

#include "MathLib/LinAlg/RowColumnIndices.h"

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
    typedef typename MathLib::RowColumnIndices<std::size_t> RowColumnIndices;
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    LocalToGlobalIndexMap(
        std::vector<LineIndex> const& rows,
        std::vector<LineIndex> const & columns)
        : _rows(rows), _columns(columns)
    {
        assert(_rows.size() == _columns.size());
        assert(!_rows.empty());
    }

    explicit LocalToGlobalIndexMap(std::vector<LineIndex> const& rows)
        : _rows(rows), _columns(rows)
    { }

    std::size_t size() const
    {
        return _rows.size();
    }

    RowColumnIndices operator[](std::size_t const mesh_item_id) const
    {
        return RowColumnIndices(_rows[mesh_item_id], _columns[mesh_item_id]);
    }

    LineIndex rowIndices(std::size_t const mesh_item_id) const
    {
        return _rows[mesh_item_id];
    }

    LineIndex columnIndices(std::size_t const mesh_item_id) const
    {
        return _columns[mesh_item_id];
    }

private:
    std::vector<LineIndex> const& _rows;
    std::vector<LineIndex> const& _columns;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALTOGLOBALINDEXMAP_H_
