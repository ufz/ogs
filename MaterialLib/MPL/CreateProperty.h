/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <boost/optional/optional.hpp>
#include <memory>
#include "PropertyType.h"

namespace BaseLib
{
class ConfigTree;
}
namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
class Property;

/// This data type is based on a std::array. It can hold pointers to objects of
/// class Property or its inheritors. The size of this array is determined by
/// the number of entries of the PropertyType enumerator.
using PropertyArray =
    std::array<std::unique_ptr<Property>, PropertyType::number_of_properties>;

/// The method reads the 'properties' tag in the prj-file and creates component
/// properties accordingly.
///
/// First, a new property iy created based on the specified property type.
/// Then, the property name is evaluated and the property is copied into the
/// properties array.
std::unique_ptr<PropertyArray> createProperties(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace MaterialPropertyLib
