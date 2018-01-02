/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
class GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface;

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class GenericNonuniformNaturalBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    template <typename Data>
    GenericNonuniformNaturalBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        unsigned const global_dim,
        std::unique_ptr<MeshLib::Mesh>&& boundary_mesh, Data&& data);

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t,
                        GlobalVector const& x,
                        GlobalMatrix& K,
                        GlobalVector& b) override;

private:
    void constructDofTable();

    /// Data used in the assembly of the specific boundary condition.
    BoundaryConditionData _data;

    std::unique_ptr<MeshLib::Mesh> _boundary_mesh;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// DOF-table (one single-component variable) of the boundary mesh.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Local assemblers for each element of the boundary mesh.
    std::vector<std::unique_ptr<
        GenericNonuniformNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;
};

}  // ProcessLib

#include "GenericNonuniformNaturalBoundaryCondition-impl.h"
