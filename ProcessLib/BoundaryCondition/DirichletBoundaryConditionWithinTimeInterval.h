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

#include "DirichletBoundaryCondition.h"

namespace NumLib
{
class TimeInterval;
}

namespace ProcessLib
{
class DirichletBoundaryConditionWithinTimeInterval final
    : public DirichletBoundaryCondition
{
public:
    DirichletBoundaryConditionWithinTimeInterval(
        double const start_time, double const end_time,
        Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    std::unique_ptr<NumLib::TimeInterval> _time_interval;
};

std::unique_ptr<DirichletBoundaryCondition>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ProcessLib::ParameterBase>>& parameters);

}  // namespace ProcessLib
