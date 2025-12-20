// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    bc_values.ids.clear();
    bc_values.values.clear();

    // convert mesh node ids to global index for the given component
    bc_values.ids.reserve(bc_values.ids.size() + _bc_values.ids.size());
    bc_values.values.reserve(bc_values.values.size() +
                             _bc_values.values.size());

    std::copy(_bc_values.ids.begin(), _bc_values.ids.end(),
              std::back_inserter(bc_values.ids));
    std::copy(_bc_values.values.begin(), _bc_values.values.end(),
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

    _bc_values.ids.clear();
    _bc_values.values.clear();

    for (auto const& l :
         MeshLib::views::meshLocations(_mesh, MeshLib::MeshItemType::Node))
    {
        const auto g_idx =
            _dof_table.getGlobalIndex(l, _variable_id, _component_id);

        if (g_idx < 0)
        {
            continue;
        }

        if ((*x[process_id])[g_idx] <= irreversibleDamage)
        {
            _bc_values.ids.emplace_back(g_idx);
            _bc_values.values.emplace_back(0.0);
        }
    }
}

void parsePhaseFieldIrreversibleDamageOracleBoundaryCondition(
    BaseLib::ConfigTree const& config)
{
    DBUG("Parsing PhaseFieldIrreversibleDamageOracleBoundaryCondition.");

    //! \ogs_file_param{prj__process_variables__process_variable__boundary_conditions__boundary_condition__type}
    config.checkConfigParameter(
        "type", "PhaseFieldIrreversibleDamageOracleBoundaryCondition");
}

std::unique_ptr<PhaseFieldIrreversibleDamageOracleBoundaryCondition>
createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id)
{
    DBUG("Constructing PhaseFieldIrreversibleDamageOracleBoundaryCondition.");

    return std::make_unique<
        PhaseFieldIrreversibleDamageOracleBoundaryCondition>(
        dof_table, mesh, variable_id, component_id);
}

}  // namespace ProcessLib
