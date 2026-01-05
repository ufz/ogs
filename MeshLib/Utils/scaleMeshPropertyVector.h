// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
