// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    triangular_matrix_type = TriangularMatrixType::Lower;
#ifdef USE_EIGEN_UNSUPPORTED
    scaling = false;
    restart = 30;
    l = 2;
    s = 4;
    angle = 0.7;
    smoothing = false;
    residualupdate = false;
#endif
}

EigenOption::SolverType EigenOption::getSolverType(
    const std::string& solver_name)
{
    if (solver_name == "CG")
    {
        return SolverType::CG;
    }
    if (solver_name == "LeastSquareCG")
    {
        return SolverType::LeastSquareCG;
    }
    if (solver_name == "BiCGSTAB")
    {
        return SolverType::BiCGSTAB;
    }
    if (solver_name == "BiCGSTABL")
    {
        return SolverType::BiCGSTABL;
    }
    if (solver_name == "IDRS")
    {
        return SolverType::IDRS;
    }
    if (solver_name == "IDRSTABL")
    {
        return SolverType::IDRSTABL;
    }
    if (solver_name == "SparseLU")
    {
        return SolverType::SparseLU;
    }
    if (solver_name == "PardisoLU")
    {
        return SolverType::PardisoLU;
    }
    if (solver_name == "GMRES")
    {
        return SolverType::GMRES;
    }

    OGS_FATAL("Unknown Eigen solver type `{:s}'", solver_name);
}

EigenOption::PreconType EigenOption::getPreconType(
    const std::string& precon_name)
{
    if (precon_name == "NONE")
    {
        return PreconType::NONE;
    }
    if (precon_name == "DIAGONAL")
    {
        return PreconType::DIAGONAL;
    }
    if (precon_name == "LeastSquareDIAGONAL")
    {
        return PreconType::LeastSquareDIAGONAL;
    }
    if (precon_name == "ILUT")
    {
        return PreconType::ILUT;
    }

    OGS_FATAL("Unknown Eigen preconditioner type `{:s}'", precon_name);
}

EigenOption::TriangularMatrixType EigenOption::getTriangularMatrixType(
    const std::string& triangular_matrix_name)
{
    if (triangular_matrix_name == "Lower")
    {
        return TriangularMatrixType::Lower;
    }
    if (triangular_matrix_name == "Upper")
    {
        return TriangularMatrixType::Upper;
    }
    if (triangular_matrix_name == "LowerUpper")
    {
        return TriangularMatrixType::LowerUpper;
    }

    OGS_FATAL("Unknown triangular matrix type `{:s}'", triangular_matrix_name);
}

std::string EigenOption::getSolverName(SolverType const solver_type)
{
    switch (solver_type)
    {
        case SolverType::CG:
            return "CG";
        case SolverType::LeastSquareCG:
            return "LeastSquareCG";
        case SolverType::BiCGSTAB:
            return "BiCGSTAB";
        case SolverType::BiCGSTABL:
            return "BiCGSTABL";
        case SolverType::IDRS:
            return "IDRS";
        case SolverType::IDRSTABL:
            return "IDRSTABL";
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
    switch (precon_type)
    {
        case PreconType::NONE:
            return "NONE";
        case PreconType::DIAGONAL:
            return "DIAGONAL";
        case PreconType::LeastSquareDIAGONAL:
            return "LeastSquareDIAGONAL";
        case PreconType::ILUT:
            return "ILUT";
    }
    return "Invalid";
}

std::string EigenOption::getTriangularMatrixName(
    TriangularMatrixType const triangular_matrix_type)
{
    switch (triangular_matrix_type)
    {
        case TriangularMatrixType::Lower:
            return "Lower";
        case TriangularMatrixType::Upper:
            return "Upper";
        case TriangularMatrixType::LowerUpper:
            return "LowerUpper";
    }
    return "Invalid";
}

}  // namespace MathLib
