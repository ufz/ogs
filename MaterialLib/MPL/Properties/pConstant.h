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
#ifndef MATERIALLIB_MPL_PROPERTIES_PCONSTANT_H_
#define MATERIALLIB_MPL_PROPERTIES_PCONSTANT_H_

#include "../mpProperty.h"

namespace MaterialPropertyLib
{
/**
 * The constant property class. This property simply retrieves the stored
 * constant value. It accepts all datatypes defined in PropertyDataType
 * (currently: double, Vector, Tensor, std::string)
 */
class Constant final : public Property
{
public:
    // Constant(){};
    explicit Constant(PropertyDataType const&);
};

}  // MaterialPropertyLib

#endif /* MATERIALLIB_MPL_PROPERTIES_PCONSTANT_H_ */
