/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

namespace MathLib
{
/// \brief Option for Eigen sparse solver
struct EigenOption final
{
    /// Solver type
    enum class SolverType : short
    {
        CG,
        LeastSquareCG,
        BiCGSTAB,
        BiCGSTABL,
        IDRS,
        IDRSTABL,
        SparseLU,
        PardisoLU,
        GMRES
    };

    /// Preconditioner type
    enum class PreconType : short
    {
        NONE,
        DIAGONAL,
        LeastSquareDIAGONAL,
        ILUT
    };

    /// triangular matrix type
    enum class TriangularMatrixType : short
    {
        Lower,
        Upper,
        LowerUpper
    };

    /// Linear solver type
    SolverType solver_type;
    /// Preconditioner type
    PreconType precon_type;
    /// Triangular Matrix Type
    TriangularMatrixType triangular_matrix_type;
    /// Maximum iteration count
    int max_iterations;
    /// Error tolerance
    double error_tolerance;
#ifdef USE_EIGEN_UNSUPPORTED
    /// Scaling the coefficient matrix and the RHS vector
    bool scaling;
    /// Restart value for the GMRES solver
    int restart;
    int l;
    int s;
    double angle;
    bool smoothing;
    bool residualupdate;
#endif

    /// Constructor
    ///
    /// Default options are CG, no preconditioner, iteration count 500 and
    /// tolerance 1e-10. Default matrix storage type is CRS.
    EigenOption();

    /// return a linear solver type from the solver name
    ///
    /// @param solver_name
    /// @return a linear solver type
    ///      If there is no solver type matched with the given name, INVALID
    ///      is returned.
    static SolverType getSolverType(const std::string& solver_name);

    /// return a preconditioner type from the name
    ///
    /// @param precon_name
    /// @return a preconditioner type
    ///      If there is no preconditioner type matched with the given name,
    ///      NONE is returned.
    static PreconType getPreconType(const std::string& precon_name);

    /// return a triangular matrix type from the name
    ///
    /// @param triangular_matrix_name
    /// @return a triangular_matrix type
    ///      If there is no triangular matrix type matched with the given name,
    ///      NONE is returned.
    static TriangularMatrixType getTriangularMatrixType(
        const std::string& triangular_matrix_name);

    /// return a linear solver name from the solver type
    static std::string getSolverName(SolverType const solver_type);

    /// return a preconditioner name from the preconditioner type
    static std::string getPreconName(PreconType const precon_type);

    /// return a triangular matrix name from the preconditioner type
    static std::string getTriangularMatrixName(
        TriangularMatrixType const triangular_matrix_type);
};

}  // namespace MathLib
