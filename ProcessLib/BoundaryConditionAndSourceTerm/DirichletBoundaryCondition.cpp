// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DirichletBoundaryCondition.h"

#include <algorithm>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "DirichletBoundaryConditionAuxiliaryFunctions.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
DirichletBoundaryCondition::DirichletBoundaryCondition(
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id)
    : _parameter(parameter),
      _bc_mesh(bc_mesh),
      _variable_id(variable_id),
      _component_id(component_id)
{
    checkParametersOfDirichletBoundaryCondition(_bc_mesh, dof_table_bulk,
                                                _variable_id, _component_id);

    std::vector<MeshLib::Node*> const& bc_nodes = _bc_mesh.getNodes();
    MeshLib::MeshSubset bc_mesh_subset(_bc_mesh, bc_nodes);

    // Create local DOF table from the BC mesh subset for the given variable
    // and component id.
    _dof_table_boundary = dof_table_bulk.deriveBoundaryConstrainedMap(
        variable_id, {component_id}, std::move(bc_mesh_subset));
}

void DirichletBoundaryCondition::getEssentialBCValues(
    const double t, GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    getEssentialBCValuesLocal(_parameter, _bc_mesh, *_dof_table_boundary,
                              _variable_id, _component_id, t, x, bc_values);
}

std::string parseDirichletBCConfig(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter("type", "Dirichlet");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__Dirichlet__parameter}
    auto const name = config.getConfigParameter<std::string>("parameter");
    DBUG("{}: parameter name {:s}", __FUNCTION__, name);
    return name;
}

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    std::string const& parameter_name, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters)
{
    DBUG("Constructing DirichletBoundaryCondition.");

    auto& parameter = ParameterLib::findParameter<double>(
        parameter_name, parameters, 1, &bc_mesh);

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

    return std::make_unique<DirichletBoundaryCondition>(
        parameter, bc_mesh, dof_table_bulk, variable_id, component_id);
}

}  // namespace ProcessLib
