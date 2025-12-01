/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
// TODO docu
/// The DirichletBoundaryCondition class describes a constant in space
/// and time Dirichlet boundary condition.
/// The expected parameter in the passed configuration is "value" which, when
/// not present defaults to zero.
class DirichletBoundaryCondition final : public BoundaryCondition
{
public:
    DirichletBoundaryCondition(
        ParameterLib::Parameter<double> const& parameter,
        MeshLib::Mesh const& bc_mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id);

    void getEssentialBCValues(
        const double t, GlobalVector const& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

private:
    ParameterLib::Parameter<double> const& _parameter;

    MeshLib::Mesh const& _bc_mesh;
    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> _dof_table_boundary;
    int const _variable_id;
    int const _component_id;
};

std::string parseDirichletBCConfig(BaseLib::ConfigTree const& config);

std::unique_ptr<DirichletBoundaryCondition> createDirichletBoundaryCondition(
    std::string const& parameter_name, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk, int const variable_id,
    int const component_id,
    const std::vector<std::unique_ptr<ParameterLib::ParameterBase>>&
        parameters);

}  // namespace ProcessLib
