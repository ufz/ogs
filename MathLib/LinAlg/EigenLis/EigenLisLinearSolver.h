/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <lis.h>

#include <vector>

#include "MathLib/LinAlg/LinearSolverOptions.h"
#include "MathLib/LinAlg/LinearSolverOptionsParser.h"

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

    bool solve(EigenMatrix &A, EigenVector& b, EigenVector &x);

private:
    bool solve(LisMatrix& A, LisVector& b, LisVector& x);
    std::string lis_options_;
};

} // MathLib
