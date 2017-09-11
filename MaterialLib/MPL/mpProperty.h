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
#ifndef MATERIALLIB_MPL_MPPROPERTY_H_
#define MATERIALLIB_MPL_MPPROPERTY_H_

#include "BaseLib/ConfigTree.h"
#include "mpEnums.h"
#include <array>
#include <boost/variant.hpp>

namespace MaterialPropertyLib
{


using Vector = std::array<double, 3>;
using Tensor = std::array<double, 9>;

using PropertyDataType = boost::variant<double, Vector, Tensor>;

class Property
{
protected:
    PropertyDataType _value;
public:
    Property();
    PropertyDataType value ();
};

using PropertyArray = std::array<Property*, number_of_property_enums>;
Property* newProperty(BaseLib::ConfigTree const&);

inline double getScalar (Property* p)
{
    return boost::get<double>(p->value());
}
template <typename T>
T getValue (T const&, Property* p)
{
    return boost::get<T>(p->value());
}

} //MaterialPropertyLib


#endif /* MATERIALLIB_MPL_MPPROPERTY_H_ */
