/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-12-05
 * \brief  Implementation of the LisOption class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LisOption.h"

namespace MathLib
{

LisOption::LisOption()
{
    solver_type = SolverType::CG;
    precon_type = PreconType::NONE;
    matrix_type = MatrixType::CRS;
    max_iterations = 500;
    error_tolerance = 1.e-10;
}

LisOption::SolverType LisOption::getSolverType(const std::string &solver_name)
{
#define RETURN_SOLVER_ENUM_IF_SAME_STRING(str, TypeName) \
    if (#TypeName==str) return SolverType::TypeName;

    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, CG);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCG);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, CGS);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCGSTAB);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCGSTABl);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, GPBiCG);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, TFQMR);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, Orthomin);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, GMRES);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, Jacobi);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, GaussSeidel);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, SOR);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCGSafe);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, CR);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCR);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, CRS);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCRSTAB);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, GPBiCR);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, BiCRSafe);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, FGMRESm);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, IDRs);
    RETURN_SOLVER_ENUM_IF_SAME_STRING(solver_name, MINRES);

    return SolverType::INVALID;
#undef RETURN_SOLVER_ENUM_IF_SAME_STRING
}

LisOption::PreconType LisOption::getPreconType(const std::string &precon_name)
{
#define RETURN_PRECOM_ENUM_IF_SAME_STRING(str, TypeName) \
    if (#TypeName==str) return PreconType::TypeName;

    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, NONE);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, JACOBI);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, ILU);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, SSOR);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, Hybrid);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, IplusS);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, SAINV);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, SAAMG);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, CroutILU);
    RETURN_PRECOM_ENUM_IF_SAME_STRING(precon_name, ILUT);

    return PreconType::NONE;
#undef RETURN_PRECOM_ENUM_IF_SAME_STRING
}

LisOption::MatrixType LisOption::getMatrixType(const std::string &matrix_name)
{
#define RETURN_MATRIX_ENUM_IF_SAME_STRING(str, TypeName) \
    if (#TypeName==str) return MatrixType::TypeName;

    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, CRS);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, CCS);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, MSR);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, DIA);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, ELL);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, JDS);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, BSR);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, BSC);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, VBR);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, COO);
    RETURN_MATRIX_ENUM_IF_SAME_STRING(matrix_name, DNS);

    return MatrixType::CRS;
#undef RETURN_MATRIX_ENUM_IF_SAME_STRING
}

} //MathLib
