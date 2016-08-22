/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERIONRESIDUAL_H
#define NUMLIB_CONVERGENCECRITERIONRESIDUAL_H

#include <boost/optional.hpp>
#include "ConvergenceCriterion.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace NumLib
{
class ConvergenceCriterionResidual final : public ConvergenceCriterion
{
public:
    explicit ConvergenceCriterionResidual(
        boost::optional<double>&& absolute_tolerance,
        boost::optional<double>&& relative_tolerance,
        MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return false; }
    bool hasResidualCheck() const override { return true; }

    void checkDeltaX(const GlobalVector& /*minus_delta_x*/,
                     GlobalVector const& /*x*/) override {}
    void checkResidual(const GlobalVector& residual) override;

    void preFirstIteration() override { _is_first_iteration = true; }
    void reset() override { _satisfied = true; _is_first_iteration = false; }
    bool isSatisfied() const { return _satisfied; }
private:
    const boost::optional<double> _abstol;
    const boost::optional<double> _reltol;
    const MathLib::VecNormType _norm_type;
    bool _satisfied = true;
    bool _is_first_iteration = true;
    double _residual_norm_0;
};

std::unique_ptr<ConvergenceCriterionResidual> createConvergenceCriterionResidual(
    BaseLib::ConfigTree const& config);

}  // NumLib

#endif  // NUMLIB_CONVERGENCECRITERIONRESIDUAL_H
