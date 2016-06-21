/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenOption.h"

namespace MathLib
{

EigenOption::EigenOption()
{
    solver_type = SolverType::SparseLU;
    precon_type = PreconType::NONE;
    max_iterations = static_cast<int>(1e6);
    error_tolerance = 1.e-16;
}

EigenOption::SolverType EigenOption::getSolverType(const std::string &solver_name)
{
#define RETURN_SOLVER_ENUM_IF_SAME_STRING(str, TypeName) \
    if (#TypeName==(str)) return SolverType::TypeName;

    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, CG);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCGSTAB);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, SparseLU);

    return SolverType::INVALID;
#undef RETURN_SOLVER_ENUM_IF_SAME_STRING
}

EigenOption::PreconType EigenOption::getPreconType(const std::string &precon_name)
{
#define RETURN_PRECOM_ENUM_IF_SAME_STRING(str, TypeName) \
    if (#TypeName==(str)) return PreconType::TypeName;

    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, NONE);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, DIAGONAL);

    return PreconType::NONE;
#undef RETURN_PRECOM_ENUM_IF_SAME_STRING
}

} //MathLib
