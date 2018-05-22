/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <map>
#include <utility>
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "Parameter.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace ProcessLib
{
template <typename T>
struct CurveScaledParameter final : public Parameter<T> {
    CurveScaledParameter(std::string const& name_,
                         MathLib::PiecewiseLinearInterpolation const& curve,
                         std::string referenced_parameter_name)
        : Parameter<T>(name_),
          _curve(curve),
          _referenced_parameter_name(std::move(referenced_parameter_name))
    {
    }

    bool isTimeDependent() const override { return true; }
    void initialize(
        std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const&
            parameters) override
    {
        _parameter =
            &findParameter<T>(_referenced_parameter_name, parameters, 0);
        _cache.resize(_parameter->getNumberOfComponents());
    }

    int getNumberOfComponents() const override
    {
        return _parameter->getNumberOfComponents();
    }

    std::vector<T> const& operator()(double const t,
                                     SpatialPosition const& pos) const override
    {
        auto const& tup = (*_parameter)(t, pos);
        auto const scaling = _curve.getValue(t);

        auto const num_comp = _parameter->getNumberOfComponents();
        for (int c = 0; c < num_comp; ++c)
        {
            _cache[c] = scaling * tup[c];
        }
        return _cache;
    }

private:
    MathLib::PiecewiseLinearInterpolation const& _curve;
    Parameter<T> const* _parameter;
    mutable std::vector<T> _cache;
    std::string const _referenced_parameter_name;
};

std::unique_ptr<ParameterBase> createCurveScaledParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // ProcessLib
