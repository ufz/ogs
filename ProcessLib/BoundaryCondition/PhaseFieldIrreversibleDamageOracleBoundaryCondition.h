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
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
class PhaseFieldIrreversibleDamageOracleBoundaryCondition final
    : public BoundaryCondition
{
public:
    PhaseFieldIrreversibleDamageOracleBoundaryCondition(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh, int const variable_id,
        int const component_id)
        : dof_table_(dof_table),
          mesh_(mesh),
          variable_id_(variable_id),
          component_id_(component_id)
    {
        if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
            component_id >=
                dof_table.getNumberOfVariableComponents(variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: ({:d}, "
                "{:d}), "
                "maximum values: ({:d}, {:d}).",
                variable_id, component_id, dof_table.getNumberOfVariables(),
                dof_table.getNumberOfVariableComponents(variable_id));
        }
    }

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, std::vector<GlobalVector*> const& x,
                     int const process_id) override;

private:
    NumLib::LocalToGlobalIndexMap const& dof_table_;
    MeshLib::Mesh const& mesh_;
    int const variable_id_;
    int const component_id_;

    NumLib::IndexValueVector<GlobalIndexType> bc_values_;
};

std::unique_ptr<PhaseFieldIrreversibleDamageOracleBoundaryCondition>
createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id);

}  // namespace ProcessLib
