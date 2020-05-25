/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    : parameter_(parameter),
      bc_mesh_(bc_mesh),
      nodes_in_bc_mesh_(bc_mesh.getNodes()),
      variable_id_(variable_id),
      component_id_(component_id),
      time_interval_(std::move(time_interval))
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
    : parameter_(parameter),
      bc_mesh_(bc_mesh),
      nodes_in_bc_mesh_(nodes_in_bc_mesh),
      variable_id_(variable_id),
      component_id_(component_id),
      time_interval_(std::move(time_interval))
{
    config(dof_table_bulk);
}

void DirichletBoundaryConditionWithinTimeInterval::config(
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk)
{
    checkParametersOfDirichletBoundaryCondition(bc_mesh_, dof_table_bulk,
                                                variable_id_, component_id_);

    std::vector<MeshLib::Node*> const& bc_nodes = bc_mesh_.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(bc_mesh_, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    dof_table_boundary_.reset(dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id_, {component_id_}, std::move(bc_mesh_subset)));
}

void DirichletBoundaryConditionWithinTimeInterval::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    if (time_interval_->contains(t))
    {
        getEssentialBCValuesLocal(parameter_, bc_mesh_, nodes_in_bc_mesh_,
                                  *dof_table_boundary_, variable_id_,
                                  component_id_, t, x, bc_values);
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
