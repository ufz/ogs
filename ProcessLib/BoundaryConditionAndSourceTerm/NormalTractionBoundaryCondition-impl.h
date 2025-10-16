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
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "NormalTractionBoundaryConditionLocalAssembler.h"
#include "ParameterLib/Utils.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Utils/CreateLocalAssemblers.h"

namespace ProcessLib
{
namespace NormalTractionBoundaryCondition
{

/// Computes the normal vector for a boundary element.
/// \param element The boundary element.
/// \param bulk_element The corresponding bulk element.
/// \tparam GlobalDim Global dimension of the problem.
/// \return The computed normal vector.
template <int GlobalDim>
Eigen::Vector3d computeElementNormal(const MeshLib::Element& element,
                                     const MeshLib::Element& bulk_element)
{
    Eigen::Vector3d normal;

    // TODO Extend to rotated 2d meshes and line elements.
    if (element.getGeomType() == MeshLib::MeshElemType::LINE)
    {
        Eigen::Vector3d const v1 = (element.getNode(1)->asEigenVector3d() -
                                    element.getNode(0)->asEigenVector3d())
                                       .normalized();
        normal[0] = -v1[1];
        normal[1] = v1[0];
        normal[2] = 0.;  // Replace the nan; only elements in
                         // xy-plane handled correctly.

        // Compute center of the bulk element to correctly orient the
        // normal.
        auto const c =
            MeshLib::getCenterOfGravity(bulk_element).asEigenVector3d();
        // Flip n if the normal points toward c:
        // <normal, c - n0> > 0.
        if (normal.dot(c - element.getNode(0)->asEigenVector3d()) > 0)
        {
            normal = -normal;
        }
    }
    else
    {
        auto const n =
            MeshLib::FaceRule::getSurfaceNormal(element).normalized();
        for (int d = 0; d < GlobalDim; ++d)
        {
            normal[d] = n[d];
        }
    }

    return normal;
}

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
    _element_normals =
        elements |
        ranges::views::transform(
            [&](const MeshLib::Element* e_ptr)
            {
                // Compute the corresponding bulk element for normal orientation
                auto const* bulk_element =
                    bulk_mesh.getElement((*bulk_element_ids)[e_ptr->getID()]);
                assert(bulk_element != nullptr);

                return computeElementNormal<GlobalDim>(*e_ptr, *bulk_element);
            }) |
        ranges::to<std::vector<Eigen::Vector3d>>();

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
