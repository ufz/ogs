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

namespace NumLib
{
std::unique_ptr<ConvergenceCriterion> createConvergenceCriterion(
    const BaseLib::ConfigTree& config)
{
    //! \ogs_file_param{process__convergence_criterion__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    if (type == "DeltaX") {
        return createConvergenceCriterionDeltaX(config);
    } else if (type == "Residual") {
        return createConvergenceCriterionResidual(config);
    } else if (type == "PerComponentDeltaX") {
        return createConvergenceCriterionPerComponentDeltaX(config);
    }

    OGS_FATAL("There is no convergence criterion of type `%s'.", type.c_str());
}

}  // NumLib
