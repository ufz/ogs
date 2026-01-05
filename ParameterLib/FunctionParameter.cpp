// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "FunctionParameter.h"

#include "BaseLib/ConfigTree.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "Function");

    std::vector<std::string> vec_expressions;

    //! \ogs_file_param{prj__parameters__parameter__Function__expression}
    for (auto const& p : config.getConfigSubtreeList("expression"))
    {
        std::string const expression_str = p.getValue<std::string>();
        vec_expressions.emplace_back(expression_str);
    }

    return std::make_unique<FunctionParameter<double>>(name, vec_expressions,
                                                       curves);
}

}  // namespace ParameterLib
