/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <boost/optional.hpp>
#include <memory>
#include <vector>

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
class Phase;
}

namespace MaterialPropertyLib
{
/// A method that parses the phase details and stores them in the private
/// _phases member.
///
/// This method creates the phases of the medium. Unlike a medium, a phase may
/// have a name. However, this is silly at the moment since this name still has
/// no effect (except of some benefits in regard of readability).
/// Phase components are required (a phase consists of at least one component).
/// Phase properties are optional. If not given, default properties are
/// assigned. These default properties average the component properties,
/// weighted by mole fraction.
std::vector<std::unique_ptr<Phase>> createPhases(
    boost::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
