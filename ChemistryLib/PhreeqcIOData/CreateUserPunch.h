/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <boost/optional/optional_fwd.hpp>
#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct UserPunch;

std::unique_ptr<UserPunch> createUserPunch(
    boost::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
