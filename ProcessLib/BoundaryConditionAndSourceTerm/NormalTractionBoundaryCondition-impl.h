/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <numeric>

#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "NormalTractionBoundaryConditionLocalAssembler.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{
template <int GlobalDim, template <typename /* shp fct */, int /* global dim */>
                         class LocalAssemblerImplementation>
NormalTractionBoundaryCondition<GlobalDim, LocalAssemblerImplementation>::
    NormalTractionBoundaryCondition(
        unsigned const integration_order, unsigned const shapefunction_order,
        MeshLib::Mesh const& bulk_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, MeshLib::Mesh const& bc_mesh,
        ParameterLib::Parameter<double> const& pressure)
    : _bc_mesh(bc_mesh),
      _integration_order(integration_order),
      _pressure(pressure)
{
    // Create component ids vector for the current variable.
    auto const& number_of_components =
        dof_table_bulk.getNumberOfVariableComponents(variable_id);
    std::vector<int> component_ids(number_of_components);
    std::iota(std::begin(component_ids), std::end(component_ids), 0);

    // BC mesh subset creation
    std::vector<MeshLib::Node*> const bc_nodes = _bc_mesh.getNodes();
    DBUG("Found {:d} nodes for Natural BCs for the variable {:d}",
         bc_nodes.size(), variable_id);

    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Compute normal vectors for each element in the boundary condition mesh.
    auto const* const bulk_element_ids = MeshLib::bulkElementIDs(_bc_mesh);
    assert(bulk_element_ids != nullptr);
    auto const& elements = _bc_mesh.getElements();
    _element_normals.resize(elements.size());
    for (std::size_t i = 0; i < elements.size(); ++i)
    {
        auto const& e = *elements[i];
        Eigen::Vector3d element_normal;

        // TODO Extend to rotated 2d meshes and line elements.
        if (e.getGeomType() == MeshLib::MeshElemType::LINE)
        {
            Eigen::Vector3d const v1 = (e.getNode(1)->asEigenVector3d() -
                                        e.getNode(0)->asEigenVector3d())
                                           .normalized();
            element_normal[0] = -v1[1];
            element_normal[1] = v1[0];
            element_normal[2] = 0.;  // Replace the nan; only elements in
                                     // xy-plane handled correctly.

            // Compute center of the bulk element to correctly orient the
            // normal.
            auto const* bulk_element =
                bulk_mesh.getElement((*bulk_element_ids)[e.getID()]);
            assert(bulk_element != nullptr);
            auto const c =
                MeshLib::getCenterOfGravity(*bulk_element).asEigenVector3d();
            // Flip n if the normal points toward c:
            // <element_normal, c - n0> > 0.
            if (element_normal.dot(c - e.getNode(0)->asEigenVector3d()) > 0)
            {
                element_normal = -element_normal;
            }
        }
        else
        {
            auto const n = MeshLib::FaceRule::getSurfaceNormal(e).normalized();
            for (int d = 0; d < GlobalDim; ++d)
            {
                element_normal[d] = n[d];
            }
        }
        _element_normals[i] = element_normal;
    }

    // Create local DOF table from the BC mesh subset for the given variable and
    // component ids.
    _dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, component_ids, std::move(bc_mesh_subset));

    BoundaryConditionAndSourceTerm::detail::createLocalAssemblers<
        GlobalDim, LocalAssemblerImplementation>(
        *_dof_table_boundary, shapefunction_order, _bc_mesh.getElements(),
        _local_assemblers, NumLib::IntegrationOrder{integration_order},
        _bc_mesh.isAxiallySymmetric(), _pressure, _element_normals);
}

template <int GlobalDim, template <typename /* shp fct */, int /* global dim */>
                         class LocalAssemblerImplementation>
void NormalTractionBoundaryCondition<GlobalDim, LocalAssemblerImplementation>::
    applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                   int const /*process_id*/, GlobalMatrix* K, GlobalVector& b,
                   GlobalMatrix* Jac)
{
    GlobalExecutor::executeMemberOnDereferenced(
        &NormalTractionBoundaryConditionLocalAssemblerInterface::assemble,
        _local_assemblers, *_dof_table_boundary, t, x, K, b, Jac);
}

template <int GlobalDim>
std::unique_ptr<NormalTractionBoundaryCondition<
    GlobalDim, NormalTractionBoundaryConditionLocalAssembler>>
createNormalTractionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    MeshLib::Mesh const& bulk_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing NormalTractionBoundaryCondition from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "NormalTraction");

    auto const parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__NormalTraction__parameter}
        config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", parameter_name);

    auto const& pressure = ParameterLib::findParameter<double>(
        parameter_name, parameters, 1, &bc_mesh);
    return std::make_unique<NormalTractionBoundaryCondition<
        GlobalDim, NormalTractionBoundaryConditionLocalAssembler>>(
        integration_order, shapefunction_order, bulk_mesh, dof_table,
        variable_id, bc_mesh, pressure);
}

}  // namespace NormalTractionBoundaryCondition
}  // namespace ProcessLib
