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

#ifndef OGS_PROPERTY_VARIABLE_TYPE_H
#define OGS_PROPERTY_VARIABLE_TYPE_H

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
    number_of_variables = 2  ///< Number of property variables.
};

const unsigned PropertyVariableNumber =
    static_cast<unsigned>(PropertyVariableType::number_of_variables);

}  // end namespace
}  // end namespace
#endif /* OGS_PROPERTY_VARIABLE_TYPE_H */
