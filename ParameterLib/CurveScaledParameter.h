/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Utils.h"

namespace ParameterLib
{
template <typename T>
struct CurveScaledParameter final : public Parameter<T>
{
    CurveScaledParameter(std::string const& name_,
                         MathLib::PiecewiseLinearInterpolation const& curve,
                         std::string referenced_parameter_name)
        : Parameter<T>(name_),
          curve_(curve),
          referenced_parameter_name_(std::move(referenced_parameter_name))
    {
    }

    bool isTimeDependent() const override { return true; }
    void initialize(
        std::vector<std::unique_ptr<ParameterBase>> const& parameters) override
    {
        parameter_ =
            &findParameter<T>(referenced_parameter_name_, parameters, 0);
        ParameterBase::mesh_ = parameter_->mesh();
    }

    int getNumberOfComponents() const override
    {
        return parameter_->getNumberOfComponents();
    }

    std::vector<T> operator()(double const t,
                              SpatialPosition const& pos) const override
    {
        // No local coordinate transformation here, which might happen twice
        // otherwise.
        assert(!this->coordinate_system_ ||
               "Coordinate system not expected to be set for curve scaled "
               "parameters.");

        auto const& tup = (*parameter_)(t, pos);
        auto const scaling = curve_.getValue(t);

        auto const num_comp = parameter_->getNumberOfComponents();
        std::vector<T> cache(num_comp);
        for (int c = 0; c < num_comp; ++c)
        {
            cache[c] = scaling * tup[c];
        }
        return cache;
    }

private:
    MathLib::PiecewiseLinearInterpolation const& curve_;
    Parameter<T> const* parameter_;
    std::string const referenced_parameter_name_;
};

std::unique_ptr<ParameterBase> createCurveScaledParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ParameterLib
