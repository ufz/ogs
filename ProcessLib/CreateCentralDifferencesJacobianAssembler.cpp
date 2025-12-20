/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-12-04 14:04:52
 */

#include "CreateCentralDifferencesJacobianAssembler.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "CentralDifferencesJacobianAssembler.h"

namespace ProcessLib
{

std::unique_ptr<AbstractJacobianAssembler>
createCentralDifferencesJacobianAssembler(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    config.checkConfigParameter("type", "CentralDifferences");

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
                    "Configuration parameter <{:s}> is deprecated and no "
                    "longer "
                    "supported.\n"
                    "(Removed in OGS version 6.5.6)\n\n"
                    "The numerical Jacobian assembler now uses absolute "
                    "perturbation values instead of relative scaling.\n\n"
                    "Migration required:\n"
                    "  Old approach:\n"
                    "    <relative_epsilons>1e-6 1e-6</relative_epsilons>\n"
                    "    <component_magnitudes>1e2 1e2</component_magnitudes>\n"
                    "  New approach:\n"
                    "    <epsilons>1e-4 1e-4</epsilons>\n"
                    " where the epsilons are equal component-wise to "
                    "relative_epsilons * component_magnitudes",
                    key);
            }
        }
    }

    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CentralDifferences__epsilons}
    auto epsilons = config.getConfigParameter<std::vector<double>>("epsilons");

    return std::make_unique<CentralDifferencesJacobianAssembler>(
        std::move(epsilons));
}

}  // namespace ProcessLib
