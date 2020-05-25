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

namespace ProcessLib
{
class GenericNaturalBoundaryConditionLocalAssemblerInterface;

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class GenericNaturalBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    template <typename Data>
    GenericNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim, MeshLib::Mesh const& bc_mesh, Data&& data);

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac) override;

private:
    /// Data used in the assembly of the specific boundary condition.
    BoundaryConditionData data_;

    /// A lower-dimensional mesh on which the boundary condition is defined.
    MeshLib::Mesh const& bc_mesh_;

    /// Local dof table, a subset of the global one restricted to the
    /// participating number of elements_ of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_boundary_;

    /// Local assemblers for each element of number of elements_.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        local_assemblers_;
};

}  // namespace ProcessLib

#include "GenericNaturalBoundaryCondition-impl.h"
