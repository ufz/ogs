/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENOPTION_H_
#define EIGENOPTION_H_

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
        BiCGSTAB,
        SparseLU,
        PardisoLU,
        GMRES
    };

    /// Preconditioner type
    enum class PreconType : short
    {
        NONE,
        DIAGONAL,
        ILUT
    };

    /// Linear solver type
    SolverType solver_type;
    /// Preconditioner type
    PreconType precon_type;
    /// Maximum iteration count
    int max_iterations;
    /// Error tolerance
    double error_tolerance;
#ifdef USE_EIGEN_UNSUPPORTED
    /// Scaling the coefficient matrix and the RHS bector
    bool scaling;
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
    static SolverType getSolverType(const std::string &solver_name);

    /// return a preconditioner type from the name
    ///
    /// @param precon_name
    /// @return a preconditioner type
    ///      If there is no preconditioner type matched with the given name, NONE
    ///      is returned.
    static PreconType getPreconType(const std::string &precon_name);

};

}
#endif // EIGENOPTION_H_
