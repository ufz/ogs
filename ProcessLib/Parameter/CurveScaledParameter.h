/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CURVESCALEDPARAMETER_H
#define PROCESSLIB_CURVESCALEDPARAMETER_H

#include <map>
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "Parameter.h"

namespace ProcessLib
{
template <typename T>
struct CurveScaledParameter final : public Parameter<T> {
    CurveScaledParameter(MathLib::PiecewiseLinearInterpolation const& curve,
                         std::string const& parameter_name)
        : _curve(curve),
          _parameter_name(parameter_name)
    {
    }

    void initialize(
        std::vector<
            std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
        override
    {
        _parameter = &findParameter<T>(_parameter_name, parameters, 0);
        _cache.resize(_parameter->getNumberOfComponents());
    }

    unsigned getNumberOfComponents() const override
    {
        return _parameter->getNumberOfComponents();
    }

    std::vector<T> const& operator()(double const t,
                                     SpatialPosition const& pos) const override
    {
        auto const& tup = (*_parameter)(t, pos);
        auto const scaling = _curve.getValue(t);

        auto const num_comp = _parameter->getNumberOfComponents();
        for (std::size_t c=0; c<num_comp; ++c) {
            _cache[c] = scaling * tup[c];
        }
        return _cache;
    }

private:
    MathLib::PiecewiseLinearInterpolation const& _curve;
    Parameter<double> const* _parameter;
    mutable std::vector<double> _cache;
    std::string const _parameter_name;
};

std::unique_ptr<ParameterBase> createCurveScaledParameter(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // ProcessLib

#endif  // PROCESSLIB_CURVESCALEDPARAMETER_H
