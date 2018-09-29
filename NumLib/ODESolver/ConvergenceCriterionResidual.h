/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional.hpp>
#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace NumLib
{
//! Convergence criterion applying a single absolute or relative tolerance to
//! the whole residual vector.
//!
//! A check of the solution increment is not done.
//! If both an absolute and a relative tolerance are specified, at least one of
//! them has to be satisfied.
class ConvergenceCriterionResidual final : public ConvergenceCriterion
{
public:
    ConvergenceCriterionResidual(boost::optional<double>&& absolute_tolerance,
                                 boost::optional<double>&& relative_tolerance,
                                 const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return true; }

    /// The function will only do diagnostic output and no actual check of the
    /// solution increment is made
    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& residual) override;

private:
    const boost::optional<double> _abstol;
    const boost::optional<double> _reltol;
    double _residual_norm_0;
};

std::unique_ptr<ConvergenceCriterionResidual> createConvergenceCriterionResidual(
    BaseLib::ConfigTree const& config);

}  // NumLib
