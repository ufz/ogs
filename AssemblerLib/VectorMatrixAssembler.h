/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
#define ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_

#include "LocalToGlobalIndexMap.h"

#include "NumLib/ODESolver/Types.h"

namespace AssemblerLib
{

// TODO doc
template<typename GlobalMatrix, typename GlobalVector,
         NumLib::ODESystemTag NLTag>
class VectorMatrixAssembler;

/// Adds result of local assembler into a global vector and a global matrix.
/// The VectorMatrixAssembler executes the local assembler for a given mesh item
/// and adds the local vector and matrix entries into the global vector and
/// the global matrix. The indices in global objects are provided by
/// the LocalToGlobalIndexMap in the construction.
template<typename GlobalMatrix, typename GlobalVector>
class VectorMatrixAssembler<GlobalMatrix, GlobalVector,
        NumLib::ODESystemTag::FirstOrderImplicitQuasilinear> final
{
public:
    VectorMatrixAssembler(
        LocalToGlobalIndexMap const& data_pos)
    : _data_pos(data_pos)
    {}

    /// Executes local assembler for the given mesh item and adds the result
    /// into the global matrix and vector.
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler_>
    void operator()(std::size_t const id,
        LocalAssembler_* const local_assembler,
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const
    {
        assert(_data_pos.size() > id);

        std::vector<GlobalIndexType> indices;

        // Local matrices and vectors will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<_data_pos.getNumComponents(); ++c)
        {
            auto const& idcs = _data_pos(id, c).rows;
            indices.reserve(indices.size() + idcs.size());
            indices.insert(indices.end(), idcs.begin(), idcs.end());
        }

        std::vector<double> local_x;
        local_x.reserve(indices.size());

        for (auto i : indices) {
            local_x.emplace_back(x.get(i));
        }

        LocalToGlobalIndexMap::RowColumnIndices const r_c_indices(
                    indices, indices);

        local_assembler->assemble(t, local_x);
        local_assembler->addToGlobal(r_c_indices, M, K, b);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};


template<typename GlobalMatrix, typename GlobalVector>
class VectorMatrixAssembler<GlobalMatrix, GlobalVector,
        NumLib::ODESystemTag::NeumannBC> final
{
public:
    VectorMatrixAssembler(
        LocalToGlobalIndexMap const& data_pos)
    : _data_pos(data_pos)
    {}

    /// Executes local assembler for the given mesh item and adds the result
    /// into the global matrix and vector.
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler_>
    void operator()(std::size_t const id,
        LocalAssembler_* const local_assembler,
        const double t, GlobalVector& b) const
    {
        // TODO I hope the changes to the VectorMatrixAssembler don't break multi-components
        assert(_data_pos.size() > id);

        // TODO Refactor: GlobalMatrix and GlobalVector are always given as
        // template params but GlobalIndexType is a global constant.
        std::vector<GlobalIndexType> indices;

        // Local matrices and vectors will always be ordered by component,
        // no matter what the order of the global matrix is.
        for (unsigned c=0; c<_data_pos.getNumComponents(); ++c)
        {
            auto const& idcs = _data_pos(id, c).rows;
            indices.reserve(indices.size() + idcs.size());
            indices.insert(indices.end(), idcs.begin(), idcs.end());
        }

        LocalToGlobalIndexMap::RowColumnIndices const r_c_indices(
                    indices, indices);

        local_assembler->assemble(t);
        local_assembler->addToGlobal(r_c_indices, b);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};


// TODO remove the dummy
template<typename GlobalMatrix, typename GlobalVector>
class VectorMatrixAssembler<GlobalMatrix, GlobalVector,
        NumLib::ODESystemTag::DUMMY> final
{
public:
    VectorMatrixAssembler(
        LocalToGlobalIndexMap const& data_pos)
    : _data_pos(data_pos)
    {}

    template <typename LocalAssembler_>
    void operator()(std::size_t const id,
        LocalAssembler_* const local_assembler) const
    {
        (void) id; (void) local_assembler;
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
