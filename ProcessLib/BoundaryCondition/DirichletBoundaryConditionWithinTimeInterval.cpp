/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DirichletBoundaryConditionWithinTimeInterval.cpp
 *
 * Created on November 26, 2018, 4:59 PM
 */
#include "DirichletBoundaryConditionWithinTimeInterval.h"

#include "DirichletBoundaryCondition.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/TimeInterval.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"

#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
DirichletBoundaryConditionWithinTimeInterval::
    DirichletBoundaryConditionWithinTimeInterval(
        std::unique_ptr<BaseLib::TimeInterval> time_interval,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _nodes_in_bc_mesh(bc_mesh.getNodes()),
      _variable_id(variable_id),
      _component_id(component_id),
      _time_interval(std::move(time_interval))
{
    config(dof_table_bulk);
}

DirichletBoundaryConditionWithinTimeInterval::
    DirichletBoundaryConditionWithinTimeInterval(
        std::unique_ptr<BaseLib::TimeInterval> time_interval,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        std::vector<MeshLib::Node*> const& nodes_in_bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _nodes_in_bc_mesh(nodes_in_bc_mesh),
      _variable_id(variable_id),
      _component_id(component_id),
      _time_interval(std::move(time_interval))
{
    config(dof_table_bulk);
}

void DirichletBoundaryConditionWithinTimeInterval::config(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(_bc_mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        _variable_id, {_component_id}, std::move(bc_mesh_subset)));
}

void DirichletBoundaryConditionWithinTimeInterval::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    if (_time_interval->contains(t))
    {
        getEssentialBCValuesLocal(_parameter, _bc_mesh, _nodes_in_bc_mesh,
                                  *_dof_table_boundary, _variable_id,
                                  _component_id, t, x, bc_values);
        return;
    }

    bc_values.ids.clear();
    bc_values.values.clear();
}

std::unique_ptr<DirichletBoundaryConditionWithinTimeInterval>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    DBUG(
        "Constructing DirichletBoundaryConditionWithinTimeInterval from "
        "config.");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "DirichletWithinTimeInterval");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s}", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &bc_mesh);

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

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__DirichletWithinTimeInterval__time_interval}
    config.peekConfigParameter<std::string>("time_interval");

    return std::make_unique<DirichletBoundaryConditionWithinTimeInterval>(
        BaseLib::createTimeInterval(config), param, bc_mesh, dof_table_bulk,
        variable_id, component_id);
}

}  // namespace ProcessLib
