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

#include <memory>
#include <vector>

#include "EigenOption.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class EigenMatrix;
class EigenVector;

class EigenLinearSolverBase;

class EigenLinearSolver final
{
public:
    /**
     * Constructor
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option      A pointer to a linear solver option. In case you omit
     *                    this argument, default settings follow those of
     *                    LisOption struct.
     */
    EigenLinearSolver(const std::string& solver_name,
                      BaseLib::ConfigTree const* const option);

    /**
     * Constructor
     * @param option Eigen linear solver options.
     */
    explicit EigenLinearSolver(std::string const& solver_name,
                               EigenOption const& option);

    ~EigenLinearSolver();

    /**
     * parse linear solvers configuration
     */
    void setOption(const BaseLib::ConfigTree& option);

    /**
     * copy linear solvers options
     */
    void setOption(const EigenOption& option) { option_ = option; }

    /**
     * get linear solver options
     */
    EigenOption& getOption() { return option_; }

    bool solve(EigenMatrix& A, EigenVector& b, EigenVector& x);

protected:
    EigenOption option_;
    std::unique_ptr<EigenLinearSolverBase> solver_;
    void setRestart();
};

}  // namespace MathLib
