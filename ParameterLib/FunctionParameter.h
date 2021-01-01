/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <exprtk.hpp>
#include <utility>
#include <vector>

#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "Parameter.h"
#include "Utils.h"

namespace ParameterLib
{
/// A parameter class evaluating functions defined by
/// user-provided mathematical expressions.
///
/// Currently, x, y, z, and t are supported as variables
/// of the functions.
template <typename T>
struct FunctionParameter final : public Parameter<T>
{
    class CurveWrapper : public exprtk::ifunction<T>
    {
    public:
        CurveWrapper(MathLib::PiecewiseLinearInterpolation const& curve)
            : exprtk::ifunction<T>(1), _curve(curve)
        {
            exprtk::disable_has_side_effects(*this);
        }
        double operator()(double const& t) override
        {
            return _curve.getValue(t);
        }

    private:
        MathLib::PiecewiseLinearInterpolation const& _curve;
    };

    using symbol_table_t = exprtk::symbol_table<T>;
    using expression_t = exprtk::expression<T>;
    using parser_t = exprtk::parser<T>;
    using error_t = exprtk::parser_error::type;

    /**
     * Constructing from a vector of expressions
     *
     * @param name        the parameter's name
     * @param vec_expression_str  a vector of mathematical expressions
     * The vector size specifies the number of components of the parameter.
     */
    FunctionParameter(
        std::string const& name,
        std::vector<std::string> const& vec_expression_str,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves)
        : Parameter<T>(name, nullptr), _vec_expression_str(vec_expression_str)
    {
        // Convert curves to function objects callable by the exprtk.
        _curves.reserve(curves.size());
        std::transform(
            begin(curves), end(curves), std::back_inserter(_curves),
            [](auto const& curve) -> std::pair<std::string, CurveWrapper> {
                return {curve.first, CurveWrapper(*curve.second)};
            });

        // Create symbol table for variables and functions.
        _symbol_table.add_constants();
        _symbol_table.create_variable("x");
        _symbol_table.create_variable("y");
        _symbol_table.create_variable("z");
        _symbol_table.create_variable("t");
        for (auto& curve : _curves)
        {
            _symbol_table.add_function(curve.first, curve.second);
        }

        // Compile expressions.
        _vec_expression.resize(_vec_expression_str.size());
        for (unsigned i = 0; i < _vec_expression_str.size(); i++)
        {
            _vec_expression[i].register_symbol_table(_symbol_table);
            parser_t parser;
            if (!parser.compile(_vec_expression_str[i], _vec_expression[i]))
            {
                OGS_FATAL("Error: {:s}\tExpression: {:s}\n",
                          parser.error(),
                          _vec_expression_str[i]);
            }
        }
    }

    bool isTimeDependent() const override { return true; }

    int getNumberOfGlobalComponents() const override
    {
        return _vec_expression.size();
    }

    std::vector<T> operator()(double const t,
                              SpatialPosition const& pos) const override
    {
        std::vector<T> cache(getNumberOfGlobalComponents());
        auto& x = _symbol_table.get_variable("x")->ref();
        auto& y = _symbol_table.get_variable("y")->ref();
        auto& z = _symbol_table.get_variable("z")->ref();
        auto& time = _symbol_table.get_variable("t")->ref();
        if (!pos.getCoordinates())
        {
            OGS_FATAL(
                "FunctionParameter: The spatial position has to be set by "
                "coordinates.");
        }
        auto const coords = pos.getCoordinates().get();
        x = coords[0];
        y = coords[1];
        z = coords[2];
        time = t;

        for (unsigned i = 0; i < _vec_expression.size(); i++)
        {
            cache[i] = _vec_expression[i].value();
        }

        if (!this->_coordinate_system)
        {
            return cache;
        }

        return this->rotateWithCoordinateSystem(cache, pos);
    }

private:
    std::vector<std::string> const _vec_expression_str;
    symbol_table_t _symbol_table;
    std::vector<expression_t> _vec_expression;
    std::vector<std::pair<std::string, CurveWrapper>> _curves;
};

std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

}  // namespace ParameterLib
