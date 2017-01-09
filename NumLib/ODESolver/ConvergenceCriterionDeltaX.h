/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERIONDELTAX_H
#define NUMLIB_CONVERGENCECRITERIONDELTAX_H

#include <boost/optional.hpp>
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
    ConvergenceCriterionDeltaX(
        boost::optional<double>&& absolute_tolerance,
        boost::optional<double>&& relative_tolerance,
        MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }

    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& /*residual*/) override {}

    void reset() override { _satisfied = true; }
    bool isSatisfied() const override { return _satisfied; }
private:
    const boost::optional<double> _abstol;
    const boost::optional<double> _reltol;
    const MathLib::VecNormType _norm_type;
    bool _satisfied = true;
};

std::unique_ptr<ConvergenceCriterionDeltaX> createConvergenceCriterionDeltaX(
    BaseLib::ConfigTree const& config);

}  // NumLib

#endif  // NUMLIB_CONVERGENCECRITERIONDELTAX_H
