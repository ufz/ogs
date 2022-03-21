/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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

    // TODO make non-optional.
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__ForwardDifferences__relative_epsilons}
    auto rel_eps = config.getConfigParameterOptional<std::vector<double>>(
        "relative_epsilons");
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__ForwardDifferences__component_magnitudes}
    auto comp_mag = config.getConfigParameterOptional<std::vector<double>>(
        "component_magnitudes");

    if (rel_eps.has_value() != comp_mag.has_value())
    {
        OGS_FATAL(
            "Either both or none of <relative_epsilons> and "
            "<component_magnitudes> have to be specified.");
    }

    std::vector<double> abs_eps;

    if (rel_eps)
    {
        if (rel_eps->size() != comp_mag->size())
        {
            OGS_FATAL(
                "The numbers of components of  <relative_epsilons> and "
                "<component_magnitudes> do not match.");
        }

        abs_eps.resize(rel_eps->size());
        for (std::size_t i = 0; i < rel_eps->size(); ++i)
        {
            abs_eps[i] = (*rel_eps)[i] * (*comp_mag)[i];
        }
    }
    else
    {
        // By default 1e-8 is used as epsilon for all components.
        // TODO: remove this default value.
        abs_eps.emplace_back(1e-8);
    }

    return std::make_unique<ForwardDifferencesJacobianAssembler>(
        std::move(abs_eps));
}

}  // namespace ProcessLib
