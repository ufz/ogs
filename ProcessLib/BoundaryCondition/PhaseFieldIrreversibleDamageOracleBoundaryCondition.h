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
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/Parameter/Parameter.h"

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
        : _dof_table(dof_table),
          _mesh(mesh),
          _variable_id(variable_id),
          _component_id(component_id)
    {
        if (variable_id >= static_cast<int>(dof_table.getNumberOfVariables()) ||
            component_id >=
                dof_table.getNumberOfVariableComponents(variable_id))
        {
            OGS_FATAL(
                "Variable id or component id too high. Actual values: (%d, "
                "%d), "
                "maximum values: (%d, %d).",
                variable_id, component_id, dof_table.getNumberOfVariables(),
                dof_table.getNumberOfVariableComponents(variable_id));
        }
    }

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:
    NumLib::LocalToGlobalIndexMap const& _dof_table;
    MeshLib::Mesh const& _mesh;
    int const _variable_id;
    int const _component_id;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;
};

std::unique_ptr<PhaseFieldIrreversibleDamageOracleBoundaryCondition>
createPhaseFieldIrreversibleDamageOracleBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh,
    int const variable_id, int const component_id);

}  // namespace ProcessLib
