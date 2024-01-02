/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on November 7, 2016, 12:14 PM
 */

#pragma once

#include <utility>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

namespace ProcessLib
{
/**
 * Fetch the nodal solutions of all coupled processes from the given vector of
 * global solutions for each process into a flat vector.
 */
std::vector<double> getCoupledLocalSolutions(
    std::vector<GlobalVector*> const& global_solutions,
    std::vector<std::vector<GlobalIndexType>> const& indices);

}  // namespace ProcessLib
