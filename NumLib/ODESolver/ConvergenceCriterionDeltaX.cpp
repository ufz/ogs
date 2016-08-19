/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionDeltaX.h"
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
ConvergenceCriterionDeltaX::ConvergenceCriterionDeltaX(
    boost::optional<double>&& absolute_tolerance,
    boost::optional<double>&& relative_tolerance,
    MathLib::VecNormType norm_type)
    : _abstol(std::move(absolute_tolerance)),
      _reltol(std::move(relative_tolerance)),
      _norm_type(norm_type)
{
    if ((!_abstol) && (!_reltol))
        OGS_FATAL(
            "At least one of absolute or relative tolerance has to be "
            "specified.");
}

void ConvergenceCriterionDeltaX::checkDeltaX(const GlobalVector& minus_delta_x,
                                             GlobalVector const& x)
{
    auto error_dx = MathLib::LinAlg::norm(minus_delta_x, _norm_type);
    auto norm_x = MathLib::LinAlg::norm(x, _norm_type);

    INFO("Convergence criterion: |dx|=%.4e, |x|=%.4e, |dx|/|x|=%.4e", error_dx,
         norm_x, error_dx / norm_x);

    bool satisfied_abs = false;
    bool satisfied_rel = false;

    if (_abstol) {
        satisfied_abs = error_dx < *_abstol;
    }
    if (_reltol) {
        satisfied_rel = error_dx < *_reltol * norm_x;
    }

    _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
}

std::unique_ptr<ConvergenceCriterionDeltaX> createConvergenceCriterionDeltaX(
    const BaseLib::ConfigTree& config)
{
    config.checkConfigParameter("type", "DeltaX");

    auto abstol = config.getConfigParameterOptional<double>("abstol");
    auto reltol = config.getConfigParameterOptional<double>("reltol");
    auto const norm_type_str =
        config.getConfigParameter<std::string>("norm_type");
    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
        OGS_FATAL("Unknown vector norm type `%s'.", norm_type_str.c_str());

    return std::unique_ptr<ConvergenceCriterionDeltaX>(
        new ConvergenceCriterionDeltaX(std::move(abstol), std::move(reltol),
                                       norm_type));
}

}  // NumLib
