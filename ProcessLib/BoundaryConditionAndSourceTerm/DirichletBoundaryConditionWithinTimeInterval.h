// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

#include "BaseLib/TimeInterval.h"
#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace BaseLib
{
class ConfigTree;
struct TimeInterval;
}  // namespace BaseLib

namespace MeshLib
{
class Node;
}

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
class DirichletBoundaryConditionWithinTimeInterval final
    : public BoundaryCondition
{
public:
    DirichletBoundaryConditionWithinTimeInterval(
        BaseLib::TimeInterval time_interval,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    void config(NumLib::LocalToGlobalIndexMap const& dof_table_bulk);

private:
    ParameterLib::Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;

    BaseLib::TimeInterval const _time_interval;
};
}  // namespace ProcessLib
