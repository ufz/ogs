// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <string>
#include <tuple>

#include "GenericNaturalBoundaryCondition.h"
#include "RobinBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct RobinBoundaryConditionConfig
{
    std::string alpha_name;
    std::string u_0_name;
    std::optional<std::string> area_parameter_name;
};

using RobinBoundaryCondition =
    GenericNaturalBoundaryCondition<RobinBoundaryConditionData,
                                    RobinBoundaryConditionLocalAssembler>;

RobinBoundaryConditionConfig parseRobinBoundaryCondition(
    BaseLib::ConfigTree const& config);

/*! Creates a new uniform Robin boundary condition from the given data.
 *
 * The Robin boundary condition is given in the form
 * \f$ \alpha \cdot [ u_0 - u(x) ] \f$,
 * where the coefficients \f$ \alpha \f$ and \f$ u_0 \f$ are obtained from the
 * \c config, and \f$ u \f$ is the unknown to which the boundary condition is
 * applied.
 *
 * The value \f$ \alpha \cdot [ u_0 - u(x) ] \f$ is a flux. It replaces the
 * integrand in the boundary integral for the variable \f$ u \f$.
 */
std::unique_ptr<RobinBoundaryCondition> createRobinBoundaryCondition(
    RobinBoundaryConditionConfig const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib
