/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_VECTORMATRIXASSEMBLER_H_
#define NUMLIB_VECTORMATRIXASSEMBLER_H_

#include "NumLib/ODESolver/Types.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace NumLib
{

namespace detail
{
inline std::vector<GlobalIndexType>
getIndices(std::size_t const id,
                    NumLib::LocalToGlobalIndexMap const& dof_table)
{
    assert(dof_table.size() > id);
    std::vector<GlobalIndexType> indices;

    // Local matrices and vectors will always be ordered by component
    // no matter what the order of the global matrix is.
    for (unsigned c = 0; c < dof_table.getNumComponents(); ++c)
    {
        auto const& idcs = dof_table(id, c).rows;
        indices.reserve(indices.size() + idcs.size());
        indices.insert(indices.end(), idcs.begin(), idcs.end());
    }

    return indices;
}

template <typename GlobalVector>
std::tuple<std::vector<double>, std::vector<GlobalIndexType>>
prepareLocalVector(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        GlobalVector const& x)
{
    auto const indices = getIndices(id, dof_table);

    std::vector<double> local_x;
    local_x.reserve(indices.size());

    for (auto i : indices)
    {
        // TODO save some function calls
        local_x.emplace_back(x[i]);
    }

    return std::make_tuple(local_x, indices);
}
}

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
         typename LocalAssembler,
         NumLib::ODESystemTag ODETag>
class VectorMatrixAssembler;


//! Specialization for first-order implicit quasi-linear systems.
template<typename GlobalMatrix, typename GlobalVector, typename LocalAssembler>
class VectorMatrixAssembler<GlobalMatrix, GlobalVector, LocalAssembler,
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
    /// \attention The index \c id is not necessarily the mesh item's id.
    void assemble(std::size_t const id,
        LocalAssembler& local_assembler,
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const
    {
        auto const x_and_indices = detail::prepareLocalVector(id, _data_pos, x);
        auto const& local_x = std::get<0>(x_and_indices);
        auto const& indices = std::get<1>(x_and_indices);
        auto const r_c_indices =
            NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

        local_assembler.assemble(t, local_x);
        local_assembler.addToGlobal(r_c_indices, M, K, b);
    }

    /// Executes assembleJacobian() of the local assembler for the
    /// given mesh item passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necessarily the mesh item's id.
    void assembleJacobian(std::size_t const id,
                          LocalAssembler& local_assembler,
                          const double t,
                          GlobalVector const& x,
                          GlobalMatrix& Jac) const
    {
        auto const x_and_indices = detail::prepareLocalVector(id, _data_pos, x);
        auto const& local_x = std::get<0>(x_and_indices);
        auto const& indices = std::get<1>(x_and_indices);
        auto const r_c_indices =
            NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

        local_assembler.assembleJacobian(t, local_x);
        local_assembler.addJacobianToGlobal(r_c_indices, Jac);
    }

    /// Executes the preTimestep() method of the local assembler for the
    /// given mesh item passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necessarily the mesh item's id.
    void preTimestep(std::size_t const id,
                     LocalAssembler& local_assembler,
                     GlobalVector const& x,
                     double const t,
                     double const delta_t) const
    {
        auto const x_and_indices = detail::prepareLocalVector(id, _data_pos, x);
        auto const& local_x = std::get<0>(x_and_indices);

        local_assembler.preTimestep(local_x, t, delta_t);
    }

    /// Executes the postTimestep() method of the local assembler for the
    /// given mesh item passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necessarily the mesh item's id.
    void postTimestep(std::size_t const id,
                      LocalAssembler& local_assembler,
                      GlobalVector const& x) const
    {
        auto const x_and_indices = detail::prepareLocalVector(id, _data_pos, x);
        auto const& local_x = std::get<0>(x_and_indices);

        local_assembler.postTimestep(local_x);
    }

    /// Executes the given callback function for the given mesh item
    /// passing the mesh item's nodal d.o.f.
    /// \attention The index \c id is not necessarily the mesh item's id.
    template <typename Callback, typename... Args>
    void passLocalVector(Callback& cb, std::size_t const id,
                         GlobalVector const& x, Args&&... args)
    {
        auto const x_and_indices = detail::prepareLocalVector(id, _data_pos, x);
        auto const& local_x = std::get<0>(x_and_indices);
        auto const& indices = std::get<1>(x_and_indices);
        auto const r_c_indices =
            NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

        cb(local_x, r_c_indices, std::forward<Args>(args)...);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};


//! Specialization used to add Neumann boundary conditions.
template<typename GlobalMatrix, typename GlobalVector, typename LocalAssembler>
class VectorMatrixAssembler<GlobalMatrix, GlobalVector, LocalAssembler,
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
    /// \attention The index \c id is not necessarily the mesh item's id.
    void assemble(std::size_t const id,
        LocalAssembler& local_assembler,
        const double t, GlobalVector& b) const
    {
        auto const indices = detail::getIndices(id, _data_pos);
        auto const r_c_indices =
            NumLib::LocalToGlobalIndexMap::RowColumnIndices(indices, indices);

        local_assembler.assemble(t);
        local_assembler.addToGlobal(r_c_indices, b);
    }

private:
    LocalToGlobalIndexMap const& _data_pos;
};

}   // namespace NumLib

#endif  // NUMLIB_VECTORMATRIXASSEMBLER_H_
