/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreateNewtonRaphsonSolverParameters.h
 *  Created on July 10, 2018, 11:32 AM
 */

#pragma once

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
struct NewtonRaphsonSolverParameters;
}

namespace MaterialLib
{
NumLib::NewtonRaphsonSolverParameters createNewtonRaphsonSolverParameters(
    BaseLib::ConfigTree const& config);
}
