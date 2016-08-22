/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenOption.h"
#include "BaseLib/Error.h"

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
    if (solver_name == "CG")
        return SolverType::CG;
    if (solver_name == "BiCGSTAB")
        return SolverType::BiCGSTAB;
    if (solver_name == "SparseLU")
        return SolverType::SparseLU;

    OGS_FATAL("Unknown Eigen solver type `%s'", solver_name.c_str());
}

EigenOption::PreconType EigenOption::getPreconType(const std::string &precon_name)
{
    if (precon_name == "NONE")
        return PreconType::NONE;
    if (precon_name == "DIAGONAL")
        return PreconType::DIAGONAL;
    if (precon_name == "ILUT")
        return PreconType::ILUT;

    OGS_FATAL("Unknown Eigen preconditioner type `%s'", precon_name.c_str());
}

} //MathLib
