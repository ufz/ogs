/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
template <int GlobalDim,
          template <typename, typename, int> class LocalAssemblerImplementation>
class NormalTractionBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    NormalTractionBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, MeshLib::Mesh const& bc_mesh,
        ParameterLib::Parameter<double> const& pressure);

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac) override;

private:
    MeshLib::Mesh const& _bc_mesh;

    /// Intersection of boundary nodes and bulk mesh subset for the
    /// variable_id/component_id pair.
    std::vector<MeshLib::Node*> _nodes_subset;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of _elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// Local assemblers for each element of number of _elements.
    std::vector<
        std::unique_ptr<NormalTractionBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    ParameterLib::Parameter<double> const& _pressure;
};

template <int GlobalDim>
std::unique_ptr<NormalTractionBoundaryCondition<
    GlobalDim, NormalTractionBoundaryConditionLocalAssembler>>
createNormalTractionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib

#include "NormalTractionBoundaryCondition-impl.h"
