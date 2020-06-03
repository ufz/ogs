/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BoundaryCondition.h"

namespace BaseLib
{
class ConfigTree;
class TimeInterval;
}

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
        std::unique_ptr<BaseLib::TimeInterval> time_interval,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    DirichletBoundaryConditionWithinTimeInterval(
        std::unique_ptr<BaseLib::TimeInterval> time_interval,
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        std::vector<MeshLib::Node*> const& nodes_in_bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    ParameterLib::Parameter<double> const& parameter_;

    MeshLib::Mesh const& bc_mesh_;
    /// Some nodes in bc_mesh_
    std::vector<MeshLib::Node*> const& nodes_in_bc_mesh_;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> dof_table_boundary_;
    int const variable_id_;
    int const component_id_;

    std::unique_ptr<BaseLib::TimeInterval const> time_interval_;

    void config(NumLib::LocalToGlobalIndexMap const& dof_table_bulk);
};

std::unique_ptr<DirichletBoundaryConditionWithinTimeInterval>
createDirichletBoundaryConditionWithinTimeInterval(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>&
        parameters);

}  // namespace ProcessLib
