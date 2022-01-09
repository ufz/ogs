/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>

#include "Parameter.h"

/// The TimeDependentHeterogeneousParameter class implements a parameter that
/// varies in time and space, i.e., it is a function \f$ (t, x) \mapsto f(t, x)
/// \in T^n \f$.
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

    /// @copydoc Parameter::getNumberOfGlobalComponents()
    int getNumberOfGlobalComponents() const override;

    bool isTimeDependent() const override;

    /// @copydoc Parameter::operator()()
    std::vector<double> operator()(double const t,
                                   SpatialPosition const& pos) const override;

    /// The TimeDependentHeterogeneousParameter depends in each time step on a
    /// parameter. Since, at construction time of the
    /// TimeDependentHeterogeneousParameter other parameter needs not to be
    /// constructed and hence can't be used to setup the object this is done
    /// later on in the initialize method.
    void initialize(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters) override;

private:
    std::vector<PairTimeParameterName> _time_parameter_name_mapping;
    std::vector<PairTimeParameter> _time_parameter_mapping;
};

std::unique_ptr<ParameterBase> createTimeDependentHeterogeneousParameter(
    std::string const& name, BaseLib::ConfigTree const& config);
}  // namespace ParameterLib
