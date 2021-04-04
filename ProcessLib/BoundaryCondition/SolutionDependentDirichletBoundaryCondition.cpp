/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SolutionDependentDirichletBoundaryCondition.h"

#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
SolutionDependentDirichletBoundaryCondition::
    SolutionDependentDirichletBoundaryCondition(
        std::string property_name,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id)
    : _bc_mesh(bc_mesh), _variable_id(variable_id), _component_id(component_id)
{
    checkParametersOfDirichletBoundaryCondition(_bc_mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset)));

    if (bc_mesh.getProperties().existsPropertyVector<double>(property_name))
    {
        OGS_FATAL(
            "Found mesh property '{:s}' in the mesh '{:s}' which is for "
            "boundary assignment. This mesh property is the built-in property "
            "of the class SolutionDependentDirichletBoundaryCondition.",
            property_name, bc_mesh.getName());
    }

    _solution_dependent_bc = MeshLib::getOrCreateMeshProperty<double>(
        const_cast<MeshLib::Mesh&>(bc_mesh), property_name,
        MeshLib::MeshItemType::Node, 1);
    _solution_dependent_bc->resize(bc_mesh.getNumberOfNodes());

    ParameterLib::SpatialPosition pos;

    auto const& nodes = bc_mesh.getNodes();
    for (std::size_t i = 0; i < _solution_dependent_bc->size(); ++i)
    {
        auto const id = nodes[i]->getID();
        pos.setNodeID(id);
        (*_solution_dependent_bc)[i] = parameter(0, pos)[0];
    }

    _parameter = std::make_unique<ParameterLib::MeshNodeParameter<double>>(
        property_name, bc_mesh, *_solution_dependent_bc);
}

void SolutionDependentDirichletBoundaryCondition::getEssentialBCValues(
    double const t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    getEssentialBCValuesLocal(*_parameter, _bc_mesh, _bc_mesh.getNodes(),
                              *_dof_table_boundary, _variable_id, _component_id,
                              t, x, bc_values);
}

void SolutionDependentDirichletBoundaryCondition::postTimestep(
    double const /*t*/, std::vector<GlobalVector*> const& x,
    int const process_id)
{
    auto const& nodes = _bc_mesh.getNodes();
    for (std::size_t i = 0; i < _solution_dependent_bc->size(); ++i)
    {
        auto const id = nodes[i]->getID();
        auto const global_index = _dof_table_boundary->getGlobalIndex(
            {_bc_mesh.getID(), MeshLib::MeshItemType::Node, id}, _variable_id,
            _component_id);

        assert(global_index >= 0);
        (*_solution_dependent_bc)[i] = x[process_id]->get(global_index);
    }
}

std::unique_ptr<SolutionDependentDirichletBoundaryCondition>
createSolutionDependentDirichletBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG(
        "Constructing SolutionDependentDirichletBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "SolutionDependentDirichlet");

    auto property_name =
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__SolutionDependentDirichlet__property_name}
        config.getConfigParameter<std::string>("property_name");

    auto& initial_value_parameter = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__SolutionDependentDirichlet__initial_value_parameter}
        config.getConfigParameter<std::string>("initial_value_parameter"),
        parameters, 1, &bc_mesh);

// In case of partitioned mesh the boundary could be empty, i.e. there is no
// boundary condition.
#ifdef USE_PETSC
    // This can be extracted to createBoundaryCondition() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createBoundaryCondition().
    if (bc_mesh.getDimension() == 0 && bc_mesh.getNumberOfNodes() == 0 &&
        bc_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    return std::make_unique<SolutionDependentDirichletBoundaryCondition>(
        std::move(property_name), initial_value_parameter, bc_mesh,
        dof_table_bulk, variable_id, component_id);
}

}  // namespace ProcessLib
