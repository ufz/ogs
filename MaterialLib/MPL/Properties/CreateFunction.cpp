// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <vector>

#include "BaseLib/ConfigTree.h"
#include "Function.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MaterialPropertyLib
{
namespace
{
/// Collect the values of a list of `<expression>` subtrees into a string
/// vector.
template <typename ConfigSubtreeList>
std::vector<std::string> parseExpressionList(
    ConfigSubtreeList const& expression_configs)
{
    std::vector<std::string> expressions;
    expressions.reserve(expression_configs.size());
    std::transform(std::begin(expression_configs), std::end(expression_configs),
                   std::back_inserter(expressions),
                   [](BaseLib::ConfigTree const& p)
                   { return p.getValue<std::string>(); });
    return expressions;
}
}  // namespace

std::unique_ptr<Function> createFunction(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Function");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create Function property {:s}.", property_name);

    //! \ogs_file_param{properties__property__Function__value}
    auto const& value_config = config.getConfigSubtree("value");

    auto const value_expressions = parseExpressionList(
        //! \ogs_file_param{properties__property__Function__value__expression}
        value_config.getConfigSubtreeList("expression"));

    // For each derivative a name of the variable and the list of expressions.
    std::vector<std::pair<std::string, std::vector<std::string>>>
        dvalue_expressions;
    //! \ogs_file_param{properties__property__Function__dvalue}
    for (auto const& dvalue_config : config.getConfigSubtreeList("dvalue"))
    {
        auto variable_name =
            //! \ogs_file_param{properties__property__Function__dvalue__variable_name}
            dvalue_config.getConfigParameter<std::string>("variable_name");

        auto expressions = parseExpressionList(
            //! \ogs_file_param{properties__property__Function__dvalue__expression}
            dvalue_config.getConfigSubtreeList("expression"));

        dvalue_expressions.emplace_back(std::move(variable_name),
                                        std::move(expressions));
    }

    // For each second derivative: two variable names and expression list.
    std::vector<D2ValueConfig> d2value_expressions;
    //! \ogs_file_param{properties__property__Function__d2value}
    for (auto const& d2value_config : config.getConfigSubtreeList("d2value"))
    {
        auto variable_names_range =
            //! \ogs_file_param{properties__property__Function__d2value__variable_name}
            d2value_config.getConfigParameterList<std::string>("variable_name");

        std::vector<std::string> variable_names(variable_names_range.begin(),
                                                variable_names_range.end());

        if (variable_names.size() != 2)
        {
            OGS_FATAL(
                "Function property '{}': each <d2value> block must contain "
                "exactly two <variable_name> entries, but {:d} were given.",
                property_name, variable_names.size());
        }

        auto expressions = parseExpressionList(
            //! \ogs_file_param{properties__property__Function__d2value__expression}
            d2value_config.getConfigSubtreeList("expression"));

        d2value_expressions.emplace_back(std::move(variable_names[0]),
                                         std::move(variable_names[1]),
                                         std::move(expressions));
    }

    return std::make_unique<MaterialPropertyLib::Function>(
        std::move(property_name),
        value_expressions,
        dvalue_expressions,
        d2value_expressions,
        curves);
}
}  // namespace MaterialPropertyLib
