/**
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "pConstant.h"

namespace MaterialPropertyLib
{
/**
 * This constructor accepts single values of any data type defined in the
 * PropertyDataType definition and sets the protected attribute _value
 * of the base class Property to that value.
 */
Constant::Constant(PropertyDataType const& v)
{
    _value = v;
};

}  // MaterialPropertyLib
