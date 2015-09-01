/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

    void setX(GLOBAL_VECTOR_* x, GLOBAL_VECTOR_* x_prev_ts)
    {
        assert(x->size() == x_prev_ts->size());
        _x = x;
        _x_prev_ts = x_prev_ts;
    }

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

        std::vector<double> localX;
        std::vector<double> localX_pts;
        LocalToGlobalIndexMap::RowColumnIndices const& indices = _data_pos[id];
        LocalToGlobalIndexMap::LineIndex remapped_rows_cols;

        const unsigned element_dof = indices.rows.size();
        auto const& x = *_x;
        auto const& x_pts = *_x_prev_ts;

        auto& mcmap = _data_pos.getMeshComponentMap();
        const unsigned num_comp = mcmap.getNumComponents();

        if (_x)         localX.reserve(element_dof);
        if (_x_prev_ts) localX_pts.reserve(element_dof);
        remapped_rows_cols.reserve(element_dof);

        // The local matrix will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<num_comp; ++c)
        {
            auto const idcs = mcmap.getIndicesForComponent(indices.rows, c);
            for (auto ip : idcs)
            {
                if (_x)         localX.emplace_back(x[ip]);
                if (_x_prev_ts) localX_pts.emplace_back(x_pts[ip]);
                remapped_rows_cols.emplace_back(ip);
            }
        }

        LocalToGlobalIndexMap::RowColumnIndices const remapped_indices(
                    remapped_rows_cols, remapped_rows_cols);

        local_assembler->assemble(localX, localX_pts);
        local_assembler->addToGlobal(_A, _rhs, remapped_indices);
    }

protected:
    GLOBAL_MATRIX_ &_A;
    GLOBAL_VECTOR_ &_rhs;
    GLOBAL_VECTOR_ *_x = nullptr;
    GLOBAL_VECTOR_ *_x_prev_ts = nullptr;
    LocalToGlobalIndexMap const& _data_pos;

    GLOBAL_VECTOR_* _secondary_variables = nullptr;
    LocalToGlobalIndexMap const* _secondary_data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
