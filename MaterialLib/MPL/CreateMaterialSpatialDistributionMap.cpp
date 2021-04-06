/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "CreateMaterialSpatialDistributionMap.h"

#include "MaterialSpatialDistributionMap.h"
#include "MeshLib/Mesh.h"

namespace MaterialPropertyLib
{
std::unique_ptr<MaterialSpatialDistributionMap>
createMaterialSpatialDistributionMap(
    std::map<int, std::shared_ptr<Medium>> const& media,
    MeshLib::Mesh const& mesh)
{
    auto const material_ids = materialIDs(mesh);

    int const max_material_id =
        !material_ids
            ? 0
            : *std::max_element(begin(*material_ids), end(*material_ids));

    if (max_material_id > static_cast<int>(media.size() - 1))
    {
        WARN(
            "The maximum value of MaterialIDs in mesh is {:d}. As the given "
            "number of porous media definitions in the project file is {:d}, "
            "the maximum value of MaterialIDs in mesh must be {:d} (index "
            "starts with zero).",
            max_material_id, media.size(), max_material_id - 1);
    }

    if (max_material_id < static_cast<int>(media.size() - 1))
    {
        WARN(
            "There are {:d} porous medium definitions in the project file but "
            "only {:d} different values in the MaterialIDs vector/data_array "
            "in the mesh.",
            media.size(), max_material_id - 1);
    }
    return std::make_unique<MaterialSpatialDistributionMap>(media,
                                                            material_ids);
}
}  // namespace MaterialPropertyLib
