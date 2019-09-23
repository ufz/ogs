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

#include <boost/optional.hpp>
#include <memory>

#include "Component.h"

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
/// The method creates components based on config subtree.
///
/// Just like a phase, a component can have a name. But, in this case, the name
/// has an important task. If a name is given, a specific component class
/// referring to that name with corresponding physical material properties is
/// created.
/// Assigning a name is optional; If no name is given, a custom component
/// without predefined properties is created.
std::vector<std::unique_ptr<Component>> createComponents(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace MaterialPropertyLib
