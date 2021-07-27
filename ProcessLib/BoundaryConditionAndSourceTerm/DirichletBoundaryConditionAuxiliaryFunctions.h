/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \brief
 *  Defines functions that are shared by DirichletBoundaryCondition
 *  and DirichletBoundaryConditionWithinTimeInterval, which avoid the way of
 *  inheritance for reducing source code duplication.
 *
 * File:   DirichletBoundaryConditionAuxiliaryFunctions.h
 *
 * Created on November 28, 2018, 11:26 AM
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
    std::vector<MeshLib::Node*> const& nodes_in_bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values);
}  // namespace ProcessLib
