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
    typename GLOBAL_VECTOR_>
class VectorMatrixAssembler
{
public:
    typedef GLOBAL_MATRIX_ GLOBAL_MATRIX;
    typedef GLOBAL_VECTOR_ GLOBAL_VECTOR;

public:
    VectorMatrixAssembler(
        GLOBAL_MATRIX_ &A,
        GLOBAL_VECTOR_ &rhs,
        LocalToGlobalIndexMap const& data_pos)
    : _A(A), _rhs(rhs), _data_pos(data_pos) {}

    ~VectorMatrixAssembler() {}

    /// Executes local assembler for the given mesh item and adds the result
    /// into the global matrix and vector.
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler_>
    void operator()(std::size_t const id,
        LocalAssembler_* const local_assembler) const
    {
        assert(_data_pos.size() > id);

        LocalToGlobalIndexMap::RowColumnIndices const& indices = _data_pos[id];

        local_assembler->assemble(indices.rows.size(), indices.columns.size());
        local_assembler->addToGlobal(_A, _rhs, indices);
    }

protected:
    GLOBAL_MATRIX_ &_A;
    GLOBAL_VECTOR_ &_rhs;
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
