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

#include "scaleMeshPropertyVector.h"

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
    std::transform(pv.begin(), pv.end(), pv.begin(),
                   [factor](auto const& v) { return v * factor; });
}
}  // namespace MeshLib
