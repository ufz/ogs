/**
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   PropertyVariableType.h
 *
 * Created on November 29, 2016, 2:28 PM
 */

#pragma once

namespace MaterialLib
{
namespace Fluid
{
/// Variable that determine the property.
enum class PropertyVariableType
{
    T = 0,                   ///< temperature.
    p = 1,                   ///< pressure.
    rho = p,                 ///< density. For some models, rho substitutes p
    C = 2,                   ///< concentration.
    number_of_variables = 3  ///< Number of property variables.
};

const unsigned PropertyVariableNumber =
    static_cast<unsigned>(PropertyVariableType::number_of_variables);

}  // end namespace
}  // end namespace
