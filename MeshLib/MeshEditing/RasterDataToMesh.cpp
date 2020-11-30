/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RasterDataToMesh.h"

#include "BaseLib/StringTools.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{

namespace RasterDataToMesh
{
static bool checkMesh(MeshLib::Mesh const& mesh)
{
    if (mesh.getDimension() > 2)
    {
        ERR("This functionality is currently only available for 2D meshes.");
        return false;
    }
    return true;
}

static double evaluatePixel(double const value, double const no_data,
                            double const replacement)
{
    if (std::abs(value - no_data) < std::numeric_limits<double>::epsilon())
    {
        return replacement;
    }
    return value;
}

bool projectToNodes(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                    double const default_replacement,
                    std::string const& array_name)
{
    if (!checkMesh(mesh))
    {
        return false;
    }

    auto& nodes = mesh.getNodes();
    auto& props = mesh.getProperties();
    std::string const name =
        BaseLib::getUniqueName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Node, 1);
    double const no_data = raster.getHeader().no_data;
    std::transform(nodes.cbegin(), nodes.cend(), std::back_inserter(*vec),
                   [&](auto const node) {
                       return evaluatePixel(raster.getValueAtPoint(*node),
                                            no_data, default_replacement);
                   });
    return true;
}

bool projectToElements(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                       double const default_replacement,
                       std::string const& array_name)
{
    if (!checkMesh(mesh))
    {
        return false;
    }

    auto& elems = mesh.getElements();
    auto& props = mesh.getProperties();
    std::string const name =
        BaseLib::getUniqueName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Cell, 1);
    double const no_data = raster.getHeader().no_data;
    std::transform(elems.cbegin(), elems.cend(), std::back_inserter(*vec),
                   [&](auto const elem) {
                       auto node = getCenterOfGravity(*elem);
                       return evaluatePixel(raster.getValueAtPoint(node),
                                            no_data, default_replacement);
                   });
    return true;
}

}  // end namespace RasterDataToMesh
}  // end namespace MeshLib
