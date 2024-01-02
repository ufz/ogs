/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <memory>
#include <vector>

#include "BoundaryCondition.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MeshLib
{
class Node;
template <typename T>
class PropertyVector;
}  // namespace MeshLib

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
struct DeactivatedSubdomainMesh;

class DeactivatedSubdomainDirichlet final : public BoundaryCondition
{
public:
    DeactivatedSubdomainDirichlet(
        MeshLib::PropertyVector<unsigned char> const& is_active,
        MathLib::PiecewiseLinearInterpolation time_interval,
        ParameterLib::Parameter<double> const& parameter,
        bool const set_outer_nodes_dirichlet_values,
        DeactivatedSubdomainMesh const& subdomain,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    void config(NumLib::LocalToGlobalIndexMap const& dof_table_bulk);

private:
    ParameterLib::Parameter<double> const& _parameter;

    DeactivatedSubdomainMesh const& _subdomain;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;

    MathLib::PiecewiseLinearInterpolation const _time_interval;
    MeshLib::PropertyVector<unsigned char> const& _is_active;

    bool const _set_outer_nodes_dirichlet_values;
};
}  // namespace ProcessLib
