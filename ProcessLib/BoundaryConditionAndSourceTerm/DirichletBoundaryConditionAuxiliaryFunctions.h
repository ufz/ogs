// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 *  Defines functions that are shared by DirichletBoundaryCondition
 *  and DirichletBoundaryConditionWithinTimeInterval, which avoid the way of
 *  inheritance for reducing source code duplication.
 */
#pragma once

#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace MeshLib
{
class Mesh;
class Node;
}  // namespace MeshLib

namespace NumLib
{
class LocalToGlobalIndexMap;
template <typename>
struct IndexValueVector;
}  // namespace NumLib

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib
{
void checkParametersOfDirichletBoundaryCondition(
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    int const variable_id,
    int const component_id);

void getEssentialBCValuesLocal(
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values);

void getEssentialBCValuesLocal(
    ParameterLib::Parameter<double> const& parameter,
    MeshLib::Mesh const& bc_mesh,
    std::vector<std::size_t> const& nodes_in_bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values);
}  // namespace ProcessLib
