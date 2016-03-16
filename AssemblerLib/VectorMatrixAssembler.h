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

namespace
{
inline AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices
getRowColumnIndices(std::size_t const id,
                    AssemblerLib::LocalToGlobalIndexMap const& dof_table,
                    std::vector<GlobalIndexType>& indices)
{
    assert(dof_table.size() > id);
    assert(indices.empty());

    // Local matrices and vectors will always be ordered by component,
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices(indices,
                                                                 indices);
}

template <typename Callback, typename GlobalVector, typename... Args>
void passLocalVector_(Callback& cb, std::size_t const id,
                      AssemblerLib::LocalToGlobalIndexMap const& dof_table,
                      GlobalVector const& x, Args&&... args)
{
    std::vector<GlobalIndexType> indices;
    auto const r_c_indices = getRowColumnIndices(id, dof_table, indices);

    std::vector<double> local_x;
    local_x.reserve(indices.size());

    for (auto i : indices)
    {
        local_x.emplace_back(x.get(i));
    }

    cb(local_x, r_c_indices, std::forward<Args>(args)...);
}
}

namespace AssemblerLib
{

/*! Calls the local assemblers of FEM processes and assembles
 *  \c GlobalMatrix'es and \c GlobalVector's.
 *
 * It optionally gets the local d.o.f. from a GlobalVector using
 * a LocalToGlobalIndexMap and passes them on to the local assembler.
 *
 * Each type of equation as flagged by the \c ODETag will have a different
 * VectorMatrixAssembler type.
 */
template<typename GlobalMatrix, typename GlobalVector,
         NumLib::ODESystemTag ODETag>
class VectorMatrixAssembler;


//! Specialization for first-order implicit quasi-linear systems.
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
    template <typename LocalAssembler>
    void operator()(std::size_t const id,
        LocalAssembler* const local_assembler,
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const
    {
        auto cb = [local_assembler](
                std::vector<double> const& local_x,
                LocalToGlobalIndexMap::RowColumnIndices const& r_c_indices,
                const double t,
                GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
        {
            local_assembler->assemble(t, local_x);
            local_assembler->addToGlobal(r_c_indices, M, K, b);
        };

        passLocalVector_(cb, id, _data_pos, x, t, M, K, b);
    }

    /// Executes the preTimestep() method of the local assembler for the
    /// given mesh item passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler>
    void preTimestep(std::size_t const id,
                     LocalAssembler* const local_assembler,
                     GlobalVector const& x) const
    {
        auto cb = [local_assembler](
            std::vector<double> const& local_x,
            LocalToGlobalIndexMap::RowColumnIndices const& /*r_c_indices*/)
        {
            local_assembler->preTimestep(local_x);
        };

        passLocalVector_(cb, id, _data_pos, x);
    }

    /// Executes the postTimestep() method of the local assembler for the
    /// given mesh item passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename LocalAssembler>
    void postTimestep(std::size_t const id,
                      LocalAssembler* const local_assembler,
                      GlobalVector const& x) const
    {
        auto cb = [local_assembler](
            std::vector<double> const& local_x,
            LocalToGlobalIndexMap::RowColumnIndices const& /*r_c_indices*/)
        {
            local_assembler->postTimestep(local_x);
        };

        passLocalVector_(cb, id, _data_pos, x);
    }

    /// Executes the given callback function for the given mesh item
    /// passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename Callback, typename... Args>
    void passLocalVector(Callback& cb, std::size_t const id,
                         GlobalVector const& x, Args&&... args)
    {
        passLocalVector_(cb, id, _data_pos, x, std::forward<Args>(args)...);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};


//! Specialization used to add Neumann boundary conditions.
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
    template <typename LocalAssembler>
    void operator()(std::size_t const id,
        LocalAssembler* const local_assembler,
        const double t, GlobalVector& b) const
    {
        std::vector<GlobalIndexType> indices;
        auto const r_c_indices = getRowColumnIndices(id, _data_pos, indices);

        local_assembler->assemble(t);
        local_assembler->addToGlobal(r_c_indices, b);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_VECTORMATRIXASSEMBLER_H_
