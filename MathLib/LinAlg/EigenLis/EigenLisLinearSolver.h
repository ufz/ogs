// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "MathLib/LinAlg/LinearSolverOptions.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"
#include "MathLib/LinAlg/Lis/LisWrapper.h"

namespace MathLib
{
class EigenVector;
class EigenMatrix;
class LisMatrix;
class LisVector;

/**
 * Linear solver using Lis library with Eigen matrix and vector objects
 */
class EigenLisLinearSolver final
{
public:
    /**
     * Constructor
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param lis_options string containing the options for LIS library; options
     * will be applied in the solve method.
     */
    EigenLisLinearSolver(std::string const& solver_name,
                         std::string const& lis_options);
    /**
     * copy linear solvers options
     * @param lis_options string containing the options for LIS library; options
     * will be applied in the solve method.
     */
    void setOption(std::string const& lis_options)
    {
        lis_options_ = lis_options;
    }

    bool solve(EigenMatrix& A, EigenVector& b, EigenVector& x);

    /// Get, if the solver can handle rectangular equation systems
    bool canSolveRectangular() const { return false; }

private:
    bool solve(LisMatrix& A, LisVector& b, LisVector& x);
    std::string lis_options_;
};

}  // namespace MathLib
