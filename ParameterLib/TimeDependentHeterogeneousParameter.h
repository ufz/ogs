/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>

#include "Parameter.h"

namespace ParameterLib
{
class TimeDependentHeterogeneousParameter final : public Parameter<double>
{
public:
    using PairTimeParameterName = std::pair<double, std::string>;
    using PairTimeParameter = std::pair<double, Parameter<double> const* const>;

    TimeDependentHeterogeneousParameter(std::string name,
                                        std::vector<PairTimeParameterName>
                                            time_parameter_name_mapping);

    int getNumberOfComponents() const override;

    bool isTimeDependent() const override;

    std::vector<double> operator()(double const t,
                                   SpatialPosition const& pos) const override;

    void initialize(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters) override;

private:
    std::vector<PairTimeParameterName> _time_parameter_name_mapping;
    std::vector<PairTimeParameter> _time_parameter_mapping;
};

std::unique_ptr<ParameterBase> createTimeDependentHeterogeneousParameter(
    std::string const& name, BaseLib::ConfigTree const& config);
}  // namespace ParameterLib
