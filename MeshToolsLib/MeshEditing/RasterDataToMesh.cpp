// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "RasterDataToMesh.h"

#include <range/v3/algorithm/transform.hpp>

#include "BaseLib/StringTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"

namespace MeshToolsLib
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

/// The point at which the raster is sampled for a mesh item.
static MathLib::Point3d const& evaluationPoint(MeshLib::Node const& node)
{
    return node;
}

static MathLib::Point3d evaluationPoint(MeshLib::Element const& element)
{
    return getCenterOfGravity(element);
}

/// Creates a new property vector of the given \c item_type on \c mesh and fills
/// it with the raster values sampled at the evaluation point of each of the
/// \c mesh_items.
template <typename MeshItems>
static bool projectToMeshItems(MeshLib::Mesh& mesh,
                               GeoLib::Raster const& raster,
                               double const default_replacement,
                               std::string const& array_name,
                               MeshLib::MeshItemType const item_type,
                               MeshItems const& mesh_items)
{
    if (!checkMesh(mesh))
    {
        return false;
    }

    auto& props = mesh.getProperties();
    std::string const name =
        BaseLib::getUniqueName(props.getPropertyVectorNames(), array_name);
    auto* const values = props.createNewPropertyVector<double>(
        name, item_type, mesh_items.size(), 1);

    double const no_data = raster.getHeader().no_data;
    ranges::transform(mesh_items, values->begin(),
                      [&](auto const item)
                      {
                          return evaluatePixel(
                              raster.getValueAtPoint(evaluationPoint(*item)),
                              no_data, default_replacement);
                      });
    return true;
}

bool projectToNodes(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                    double const default_replacement,
                    std::string const& array_name)
{
    return projectToMeshItems(mesh, raster, default_replacement, array_name,
                              MeshLib::MeshItemType::Node, mesh.getNodes());
}

bool projectToElements(MeshLib::Mesh& mesh, GeoLib::Raster const& raster,
                       double const default_replacement,
                       std::string const& array_name)
{
    return projectToMeshItems(mesh, raster, default_replacement, array_name,
                              MeshLib::MeshItemType::Cell, mesh.getElements());
}

}  // end namespace RasterDataToMesh
}  // namespace MeshToolsLib
