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

#include <utility>
#include <vector>

#include <exprtk.hpp>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

#include "Parameter.h"
#include "Utils.h"

namespace ParameterLib
{
/// A parameter class evaluating functions defined by
/// user-provided mathematical expressions.
///
/// Currently, x, y, and z are supported as variables
/// of the functions.
template <typename T>
struct FunctionParameter final : public Parameter<T>
{
    using symbol_table_t = exprtk::symbol_table<T>;
    using expression_t = exprtk::expression<T>;
    using parser_t = exprtk::parser<T>;
    using error_t = exprtk::parser_error::type;

    /**
     * Constructing from a vector of expressions
     *
     * @param name        the parameter's name
     * @param mesh        the parameter's domain of definition.
     * @param vec_expression_str  a vector of mathematical expressions
     * The vector size specifies the number of components of the parameter.
     */
    FunctionParameter(std::string const& name,
                      MeshLib::Mesh const& mesh,
                      std::vector<std::string> const& vec_expression_str)
        : Parameter<T>(name, &mesh), vec_expression_str_(vec_expression_str)
    {
        symbol_table_.add_constants();
        symbol_table_.create_variable("x");
        symbol_table_.create_variable("y");
        symbol_table_.create_variable("z");

        vec_expression_.resize(vec_expression_str_.size());
        for (unsigned i = 0; i < vec_expression_str_.size(); i++)
        {
            vec_expression_[i].register_symbol_table(symbol_table_);
            parser_t parser;
            if (!parser.compile(vec_expression_str_[i], vec_expression_[i]))
            {
                OGS_FATAL("Error: {:s}\tExpression: {:s}\n",
                          parser.error(),
                          vec_expression_str_[i]);
            }
        }
    }

    bool isTimeDependent() const override { return false; }

    int getNumberOfComponents() const override
    {
        return vec_expression_.size();
    }

    std::vector<T> operator()(double const /*t*/,
                              SpatialPosition const& pos) const override
    {
        std::vector<T> cache(getNumberOfComponents());
        auto& x = symbol_table_.get_variable("x")->ref();
        auto& y = symbol_table_.get_variable("y")->ref();
        auto& z = symbol_table_.get_variable("z")->ref();
        if (pos.getCoordinates())
        {
            auto const coords = pos.getCoordinates().get();
            x = coords[0];
            y = coords[1];
            z = coords[2];
        }
        else if (pos.getNodeID())
        {
            auto const& node =
                *ParameterBase::mesh_->getNode(pos.getNodeID().get());
            x = node[0];
            y = node[1];
            z = node[2];
        }

        for (unsigned i = 0; i < vec_expression_.size(); i++)
        {
            cache[i] = vec_expression_[i].value();
        }

        if (!this->coordinate_system_)
        {
            return cache;
        }

        return this->rotateWithCoordinateSystem(cache, pos);
    }

private:
    std::vector<std::string> const vec_expression_str_;
    symbol_table_t symbol_table_;
    std::vector<expression_t> vec_expression_;
};

std::unique_ptr<ParameterBase> createFunctionParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh);

}  // namespace ParameterLib
