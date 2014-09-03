/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
#define ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_

#include "LocalToGlobalIndexMap.h"

namespace AssemblerLib
{

/// Adds result of local assembler into a global vector and a global matrix.
/// The VectorMatrixAssembler executes the local assembler for a given mesh item
/// and adds the local vector and matrix entries into the global vector and
/// the global matrix. The indices in global objects are provided by
/// the LocalToGlobalIndexMap in the construction.
template<
    typename GLOBAL_MATRIX_,
    typename GLOBAL_VECTOR_,
    typename MESH_ITEM_,
    typename ASSEMBLER_,
    typename LOCAL_MATRIX_,
    typename LOCAL_VECTOR_>
class VectorMatrixAssembler
{
public:
    typedef GLOBAL_MATRIX_ GLOBAL_MATRIX;
    typedef GLOBAL_VECTOR_ GLOBAL_VECTOR;
    typedef LOCAL_MATRIX_ LOCAL_MATRIX;
    typedef LOCAL_VECTOR_ LOCAL_VECTOR;

public:
    VectorMatrixAssembler(
        GLOBAL_MATRIX_ &A,
        GLOBAL_VECTOR_ &rhs,
        ASSEMBLER_ &local_assembler,
        LocalToGlobalIndexMap const& data_pos)
    : _A(A), _rhs(rhs), _local_assembler(local_assembler), _data_pos(data_pos) {}

    ~VectorMatrixAssembler() {}

    /// Executes local assembler for the given mesh item and adds the result
    /// into the global matrix and vector.
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    void operator()(const MESH_ITEM_* item, std::size_t id) const
    {
        assert(_data_pos.size() > id);

        LocalToGlobalIndexMap::RowColumnIndices const& indices = _data_pos[id];
        LOCAL_MATRIX_ local_A(indices.rows.size(), indices.columns.size());
        LOCAL_VECTOR_ local_rhs(indices.rows.size());

        _local_assembler(*item, local_A, local_rhs);
        _A.add(indices, local_A);
        _rhs.add(indices.rows, local_rhs);
    }

protected:
    GLOBAL_MATRIX_ &_A;
    GLOBAL_VECTOR_ &_rhs;
    ASSEMBLER_ &_local_assembler;
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
