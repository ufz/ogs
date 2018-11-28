/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace MeshLib
{
class Mesh;
}

namespace NumLib
{
class LocalToGlobalIndexMap;
template <typename>
struct IndexValueVector;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

void checkParametersOfDirichletBoundaryCondition(
    MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    int const variable_id,
    int const component_id);

void getEssentialBCValuesLocal(
    Parameter<double> const& parameter, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
    int const variable_id, int const component_id, const double t,
    GlobalVector const& x,
    NumLib::IndexValueVector<GlobalIndexType>& bc_values);
}  // end of name space
