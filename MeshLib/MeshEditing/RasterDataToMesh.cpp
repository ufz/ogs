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

static std::string getValidName(std::vector<std::string> const& vec_names,
                                std::string const& array_name)
{
    std::string new_name = array_name;
    std::size_t count = 1;
    while (std::find(vec_names.cbegin(), vec_names.cend(), new_name) !=
           vec_names.end())
    {
        count++;
        new_name = array_name + "-" + std::to_string(count);
    }
    return new_name;
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
                    double const def, std::string const& array_name)
{
    if (!checkMesh(mesh))
    {
        return false;
    }

    auto& nodes = mesh.getNodes();
    auto& props = mesh.getProperties();
    std::string const name =
        getValidName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Node, 1);
    double const no_data = raster.getHeader().no_data;
    std::transform(nodes.cbegin(), nodes.cend(), std::back_inserter(*vec), [&](auto const node)
    {
        return evaluatePixel(raster.getValueAtPoint(*node), no_data, def);
    });
    return true;
}

bool projectToElements(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                       double const def, std::string const& array_name)
{
    if (!checkMesh(mesh))
    {
        return false;
    }

    auto& elems = mesh.getElements();
    auto& props = mesh.getProperties();
    std::string const name =
        getValidName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Cell, 1);
    double const no_data = raster.getHeader().no_data;
    std::transform(elems.cbegin(), elems.cend(), std::back_inserter(*vec), [&](auto const elem) {
        auto node = elem->getCenterOfGravity();
        return evaluatePixel(raster.getValueAtPoint(node), no_data, def);
    });
    return true;
}

}  // end namespace RasterDataToMesh
}  // end namespace MeshLib
