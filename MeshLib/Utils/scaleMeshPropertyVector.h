/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>

namespace MeshLib
{
class Mesh;

/// Scales the mesh property with name \c property_name by given \c factor.
/// \note The property must be a "double" property.
void scaleMeshPropertyVector(Mesh& mesh,
                             std::string const& property_name,
                             double factor);
}  // namespace MeshLib
