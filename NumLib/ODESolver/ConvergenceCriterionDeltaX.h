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

#include <optional>

#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace NumLib
{
//! Convergence criterion applying a single absolute or relative tolerance to
//! the whole solution increment vector.
//!
//! A residual check is not done.
//! If both an absolute and a relative tolerance are specified, at least one of
//! them has to be satisfied.
class ConvergenceCriterionDeltaX final : public ConvergenceCriterion
{
public:
    ConvergenceCriterionDeltaX(std::optional<double>&& absolute_tolerance,
                               std::optional<double>&& relative_tolerance,
                               const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }
    bool hasNonNegativeDamping() const override { return false; }

    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& /*residual*/) override {}
    double getDampingFactor(GlobalVector const& /*minus_delta_x*/,
                            GlobalVector const& /*x*/,
                            double damping_scalar) override
    {
        return damping_scalar;
    }

    void reset() override { this->_satisfied = true; }

private:
    const std::optional<double> _abstol;
    const std::optional<double> _reltol;
};

std::unique_ptr<ConvergenceCriterionDeltaX> createConvergenceCriterionDeltaX(
    BaseLib::ConfigTree const& config);

}  // namespace NumLib
