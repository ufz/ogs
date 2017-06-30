/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    const MathLib::VecNormType norm_type)
    : ConvergenceCriterion(norm_type), _abstol(std::move(absolute_tolerance)),
      _reltol(std::move(relative_tolerance))
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
        satisfied_rel = checkRelativeTolerance(*_reltol, error_dx, norm_x);
    }

    _satisfied = _satisfied && (satisfied_abs || satisfied_rel);
}

std::unique_ptr<ConvergenceCriterionDeltaX> createConvergenceCriterionDeltaX(
    const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__type}
    config.checkConfigParameter("type", "DeltaX");

    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__DeltaX__abstol}
    auto abstol = config.getConfigParameterOptional<double>("abstol");
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__DeltaX__reltol}
    auto reltol = config.getConfigParameterOptional<double>("reltol");
    auto const norm_type_str =
        //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__DeltaX__norm_type}
        config.getConfigParameter<std::string>("norm_type");
    auto const norm_type = MathLib::convertStringToVecNormType(norm_type_str);

    if (norm_type == MathLib::VecNormType::INVALID)
        OGS_FATAL("Unknown vector norm type `%s'.", norm_type_str.c_str());

    return std::make_unique<ConvergenceCriterionDeltaX>(
        std::move(abstol), std::move(reltol), norm_type);
}

}  // NumLib
