// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "scaleMeshPropertyVector.h"

#include <range/v3/algorithm/transform.hpp>

#include "MeshLib/Mesh.h"

namespace MeshLib
{
void scaleMeshPropertyVector(MeshLib::Mesh& mesh,
                             std::string const& property_name,
                             double factor)
{
    if (!mesh.getProperties().existsPropertyVector<double>(property_name))
    {
        WARN("Did not find PropertyVector '{:s}' for scaling.", property_name);
        return;
    }
    auto& pv = *mesh.getProperties().getPropertyVector<double>(property_name);
    ranges::transform(pv, pv.begin(),
                      [factor](auto const& v) { return v * factor; });
}
}  // namespace MeshLib
