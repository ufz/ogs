/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   DirichletBoundaryConditionWithinTimeInterval.h
 *
 * Created on November 26, 2018, 4:59 PM
 */
#pragma once

#include <memory>

#include "BoundaryCondition.h"

namespace BaseLib
{
class ConfigTree;
class TimeInterval;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

class DirichletBoundaryConditionWithinTimeInterval final
    : public BoundaryCondition
{
public:
    DirichletBoundaryConditionWithinTimeInterval(
        std::unique_ptr<BaseLib::TimeInterval> time_interval,
        Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;

    std::unique_ptr<BaseLib::TimeInterval> _time_interval;
};

std::unique_ptr<DirichletBoundaryConditionWithinTimeInterval>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
