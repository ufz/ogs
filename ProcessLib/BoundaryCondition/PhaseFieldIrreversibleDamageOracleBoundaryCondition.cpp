/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PhaseFieldIrreversibleDamageOracleBoundaryCondition.h"

#include <algorithm>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"

namespace ProcessLib
{
void PhaseFieldIrreversibleDamageOracleBoundaryCondition::getEssentialBCValues(
    const double /*t*/, GlobalVector const& /*x*/,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values) const
{
    ParameterLib::SpatialPosition pos;

    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + bc_values_.ids.size());
    bc_values.values.reserve(bc_values.values.size() +
                             bc_values_.values.size());

    std::copy(bc_values_.ids.begin(), bc_values_.ids.end(),
              std::back_inserter(bc_values.ids));
    std::copy(bc_values_.values.begin(), bc_values_.values.end(),
              std::back_inserter(bc_values.values));
}

// update new values and corresponding indices.
void PhaseFieldIrreversibleDamageOracleBoundaryCondition::preTimestep(
    const double /*t*/, std::vector<GlobalVector*> const& x,
    int const process_id)
{
    // phase-field variable is considered irreversible if it loses more than 95%
    // of the stiffness, which is a widely used threshold.
    double irreversibleDamage = 0.05;

    bc_values_.ids.clear();
    bc_values_.values.clear();

    auto const mesh_id = mesh_.getID();
    auto const& nodes = mesh_.getNodes();
    for (auto const* n : nodes)
    {
        std::size_t node_id = n->getID();
        MeshLib::Location l(mesh_id, MeshLib::MeshItemType::Node, node_id);
        const auto g_idx =
            dof_table_.getGlobalIndex(l, variable_id_, component_id_);

        if ((*x[process_id])[node_id] <= irreversibleDamage)
        {
            bc_values_.ids.emplace_back(g_idx);
            bc_values_.values.emplace_back(0.0);
        }
    }
}

std::unique_ptr<PhaseFieldIrreversibleDamageOracleBoundaryCondition>
createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id)
{
    DBUG(
        "Constructing PhaseFieldIrreversibleDamageOracleBoundaryCondition from "
        "config.");
    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter(
        "type", "PhaseFieldIrreversibleDamageOracleBoundaryCondition");

    return std::make_unique<
        PhaseFieldIrreversibleDamageOracleBoundaryCondition>(
        dof_table, mesh, variable_id, component_id);
}

}  // namespace ProcessLib
