/**
 * \file
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/MeshSubset.h"
#include "BoundaryCondition.h"
#include "PressureBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
namespace PressureBoundaryCondition
{
class PressureBoundaryConditionLocalAssemblerInterface;

template <template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class PressureBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    PressureBoundaryCondition(
        bool const is_axially_symmetric, unsigned const integration_order,
        unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, unsigned const global_dim,
        std::vector<MeshLib::Element*>&& elements,
        Parameter<double> const& pressure);

    ~PressureBoundaryCondition() override;

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t,
                        GlobalVector const& x,
                        GlobalMatrix& K,
                        GlobalVector& b) override;

private:
    /// Vector of lower-dimensional elements on which the boundary condition is
    /// defined.
    std::vector<MeshLib::Element*> _elements;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<
        std::unique_ptr<PressureBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    Parameter<double> const& _pressure;
};

std::unique_ptr<
    PressureBoundaryCondition<PressureBoundaryConditionLocalAssembler>>
createPressureBoundaryCondition(
    BaseLib::ConfigTree const& config,
    std::vector<MeshLib::Element*>&& elements,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    bool is_axially_symmetric, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

}  // namespace PressureBoundaryCondition
}  // namespace ProcessLib

#include "PressureBoundaryCondition-impl.h"
