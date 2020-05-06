/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Constant.h"

namespace MaterialPropertyLib
{
Constant::Constant(PropertyDataType const& v)
{
    value_ = v;
    dvalue_ = std::visit(
        [](auto const& value) -> PropertyDataType { return decltype(value){}; },
        v);
};
}  // namespace MaterialPropertyLib
