/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateForwardDifferencesJacobianAssembler.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ForwardDifferencesJacobianAssembler.h"

namespace ProcessLib
{
std::unique_ptr<AbstractJacobianAssembler>
createForwardDifferencesJacobianAssembler(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    config.checkConfigParameter("type", "ForwardDifferences");

    // TODO: to be removed
    {
        std::vector<std::string> const deprecated_keys{"relative_epsilons",
                                                       "component_magnitudes"};
        for (auto const& key : deprecated_keys)
        {
            auto const deprecated_parameter =
                //! \ogs_file_special
                config.getConfigParameterOptional<std::vector<double>>(key);

            if ((deprecated_parameter.has_value() &&
                 !deprecated_parameter->empty()))
            {
                OGS_FATAL(
                    "Using {:s} is not allowed anymore as of r6.5.6. Please "
                    "replace <relative_epsilons> and <component_magnitudes> "
                    "with <epsilons> to specify variable component "
                    "perturbations (excluding deformation).",
                    key);
            }
        }
    }

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__ForwardDifferences__epsilons}
    auto epsilons = config.getConfigParameter<std::vector<double>>("epsilons");

    return std::make_unique<ForwardDifferencesJacobianAssembler>(
        std::move(epsilons));
}

}  // namespace ProcessLib
