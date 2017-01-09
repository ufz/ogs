/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
#ifdef USE_EIGEN_UNSUPPORTED
    scaling = false;
#endif
}

EigenOption::SolverType EigenOption::getSolverType(const std::string &solver_name)
{
    if (solver_name == "CG")
        return SolverType::CG;
    if (solver_name == "BiCGSTAB")
        return SolverType::BiCGSTAB;
    if (solver_name == "SparseLU")
        return SolverType::SparseLU;
    if (solver_name == "PardisoLU")
        return SolverType::PardisoLU;
    if (solver_name == "GMRES")
        return SolverType::GMRES;

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

std::string EigenOption::getSolverName(SolverType const solver_type)
{
    switch (solver_type) {
        case SolverType::CG:
            return "CG";
        case SolverType::BiCGSTAB:
            return "BiCGSTAB";
        case SolverType::SparseLU:
            return "SparseLU";
        case SolverType::PardisoLU:
            return "PardisoLU";
        case SolverType::GMRES:
            return "GMRES";
    }
    return "Invalid";
}

std::string EigenOption::getPreconName(PreconType const precon_type)
{
    switch (precon_type) {
        case PreconType::NONE:
            return "NONE";
        case PreconType::DIAGONAL:
            return "DIAGONAL";
        case PreconType::ILUT:
            return "ILUT";
    }
    return "Invalid";
}

} //MathLib
