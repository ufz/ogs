/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterionResidual.h"
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace NumLib
{
ConvergenceCriterionResidual::ConvergenceCriterionResidual(
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

void ConvergenceCriterionResidual::checkResidual(const GlobalVector& residual)
{
    auto norm_res = MathLib::LinAlg::norm(residual, _norm_type);

    if (_is_first_iteration) {
        INFO("Convergence criterion: |r0|=%.4e", norm_res);
        _residual_norm_0 = norm_res;
    } else {
        INFO("Convergence criterion: |r|=%.4e |r0|=%.4e |r|/|r0|=%.4e",
             norm_res, _residual_norm_0, norm_res / _residual_norm_0);
    }

    bool satisfied_abs = false;
    bool satisfied_rel = false;

    if (_abstol) {
        satisfied_abs = norm_res < *_abstol;
    }
    if (_reltol && !_is_first_iteration) {
        satisfied_rel =
            checkRelativeTolerance(*_reltol, norm_res, _residual_norm_0);
    }

    _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
}

std::unique_ptr<ConvergenceCriterionResidual>
createConvergenceCriterionResidual(const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{process__convergence_criterion__type}
    config.checkConfigParameter("type", "Residual");

    //! \ogs_file_param{process__convergence_criterion__Residual__abstol}
    auto abstol = config.getConfigParameterOptional<double>("abstol");
    //! \ogs_file_param{process__convergence_criterion__Residual__reltol}
    auto reltol = config.getConfigParameterOptional<double>("reltol");
    auto const norm_type_str =
        //! \ogs_file_param{process__convergence_criterion__Residual__norm_type}
        config.getConfigParameter<std::string>("norm_type");
    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
        OGS_FATAL("Unknown vector norm type `%s'.", norm_type_str.c_str());

    return std::unique_ptr<ConvergenceCriterionResidual>(
        new ConvergenceCriterionResidual(std::move(abstol), std::move(reltol),
                                         norm_type));
}

}  // NumLib
