/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BoundaryCondition.h"
#include "MeshLib/MeshSubset.h"
#include "NormalTractionBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{
class NormalTractionBoundaryConditionLocalAssemblerInterface;

/// The normal traction boundary condition is a special type of Neumann boundary
/// condition where the given value is applied in the direction of the element's
/// normal vector \f$\mathbf{n}\f$:
/// \f[
///      \bar{t} := \sigma \mathbf{n} = p \mathbf{n},
/// \f]
/// where \f$p\f$ is the value on the boundary given by the parameter tag.
template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class NormalTractionBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    NormalTractionBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, unsigned const global_dim,
        MeshLib::Mesh const& bc_mesh,
        ParameterLib::Parameter<double> const& pressure);

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac) override;

private:
    MeshLib::Mesh const& bc_mesh_;

    /// Intersection of boundary nodes and bulk mesh subset for the
    /// variable_id/component_id pair.
    std::vector<MeshLib::Node*> nodes_subset_;

    std::unique_ptr<MeshLib::MeshSubset const> mesh_subset_all_nodes_;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements_ of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_boundary_;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const integration_order_;

    /// Local assemblers for each element of number of elements_.
    std::vector<
        std::unique_ptr<NormalTractionBoundaryConditionLocalAssemblerInterface>>
        local_assemblers_;

    ParameterLib::Parameter<double> const& pressure_;
};

std::unique_ptr<NormalTractionBoundaryCondition<
    NormalTractionBoundaryConditionLocalAssembler>>
createNormalTractionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib

#include "NormalTractionBoundaryCondition-impl.h"
