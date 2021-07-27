/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
#include <vector>

#include "BaseLib/TimeInterval.h"
#include "BoundaryCondition.h"

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
