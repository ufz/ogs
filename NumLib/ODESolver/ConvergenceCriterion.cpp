/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConvergenceCriterion.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ConvergenceCriterionDeltaX.h"
#include "ConvergenceCriterionResidual.h"
#include "ConvergenceCriterionPerComponentDeltaX.h"
#include "ConvergenceCriterionPerComponentResidual.h"

namespace NumLib
{
std::unique_ptr<ConvergenceCriterion> createConvergenceCriterion(
    const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__convergence_criterion__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "DeltaX") {
        return createConvergenceCriterionDeltaX(config);
    } else if (type == "Residual") {
        return createConvergenceCriterionResidual(config);
    } else if (type == "PerComponentDeltaX") {
        return createConvergenceCriterionPerComponentDeltaX(config);
    } else if (type == "PerComponentResidual") {
        return createConvergenceCriterionPerComponentResidual(config);
    }

    OGS_FATAL("There is no convergence criterion of type `%s'.", type.c_str());
}

bool checkRelativeTolerance(const double reltol, const double numerator,
                            const double denominator)
{
    auto const eps = std::numeric_limits<double>::epsilon();
    return std::abs(numerator) <
           std::abs(reltol) * (std::abs(denominator) + eps);
}

}  // NumLib
