// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "BaseLib/OgsAsmThreads.h"
#include "FunctionEvaluation.h"
#include "Parameter.h"

namespace ParameterLib
{
/// A parameter class evaluating functions defined by
/// user-provided mathematical expressions.
///
/// Currently, x, y, z, and t are supported as variables
/// of the functions.
///
/// Uses per-thread storage for thread-safe evaluation.
template <typename T>
struct FunctionParameter final : public Parameter<T>
{
    /**
     * Constructing from a vector of expressions
     *
     * \param name        the parameter's name
     * \param vec_expression_str  a vector of mathematical expressions
     * \param curves      named list of curves used by expressions.
     * The vector size specifies the number of components of the parameter.
     */
    FunctionParameter(
        std::string const& name,
        std::vector<std::string> const& vec_expression_str,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves)
        : Parameter<T>(name, nullptr),
          function_evaluation_(BaseLib::getNumberOfThreads(),
                               {"x", "y", "z", "t"}, vec_expression_str, curves)
    {
    }

    bool isTimeDependent() const override
    {
        return function_evaluation_.isTimeDependent();
    }

    int getNumberOfGlobalComponents() const override
    {
        return function_evaluation_.getNumberOfComponents();
    }

    std::vector<T> operator()(double const t,
                              SpatialPosition const& pos) const override
    {
        std::vector<T> values(static_cast<std::size_t>(
            function_evaluation_.getNumberOfComponents()));
        function_evaluation_.evaluate(pos, t, values);
        if (!this->_coordinate_system)
        {
            return values;
        }
        return this->rotateWithCoordinateSystem(values, pos);
    }

    FunctionEvaluation function_evaluation_;
};

std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ParameterLib
