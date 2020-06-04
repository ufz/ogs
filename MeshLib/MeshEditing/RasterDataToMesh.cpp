/**
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

bool checkMesh2D(MeshLib::Mesh const& mesh)
{
    if (mesh.getDimension() > 2)
    {
        ERR("This functionality is currently only available for 2D meshes.");
        return false;
    }
    return true;
}
std::string getValidName(std::vector<std::string> const& vec_names,
                     std::string const& array_name)
{
    std::string new_name = array_name;
    std::size_t count = 1;
    bool name_exists = std::find(vec_names.cbegin(), vec_names.cend(),
                                 array_name) != vec_names.end();
    while (name_exists)
    {
        count++;
        new_name = array_name + std::to_string(count);
        name_exists = std::find(vec_names.cbegin(), vec_names.cend(),
                                array_name) != vec_names.end();
    }
    return new_name;
}

double evaluatePixel(double value, double no_data, double replacement)
{
    if (std::abs(value - no_data) < std::numeric_limits<double>::epsilon())
        return replacement;
    return value;
}

bool projectToNodes(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                    double const default, std::string const& array_name)
{
    if (!checkMesh2D(mesh))
        return false;

    auto& nodes = mesh.getNodes();
    auto& props = mesh.getProperties();
    std::string const name =
        getValidName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Node, 1);
    double const no_data = raster.getHeader().no_data;
    for (auto node : nodes)
    {
        vec->push_back(
            evaluatePixel(raster.getValueAtPoint(*node), no_data, default));
    }
    return true;
}

bool projectToElements(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                       double const default, std::string const& array_name)
{
    if (!checkMesh2D(mesh))
        return false;

    auto& elems = mesh.getElements();
    auto& props = mesh.getProperties();
    std::string const name =
        getValidName(props.getPropertyVectorNames(), array_name);
    auto vec = props.createNewPropertyVector<double>(
        name, MeshLib::MeshItemType::Cell, 1);
    double const no_data = raster.getHeader().no_data;
    for (auto elem : elems)
    {
        auto node = elem->getCenterOfGravity();
        vec->push_back(
            evaluatePixel(raster.getValueAtPoint(node), no_data, default));
    }
    return true;
}

}  // namespace RasterDataToMesh
} // end namespace MeshLib
