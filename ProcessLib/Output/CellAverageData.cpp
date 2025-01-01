/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CellAverageData.h"

#include "MeshLib/Utils/getOrCreateMeshProperty.h"

namespace ProcessLib
{
MeshLib::PropertyVector<double>& CellAverageData::getOrCreatePropertyVector(
    std::string const& name, unsigned const num_comp)
{
    if (auto const it = cell_averages_.find(name); it != cell_averages_.end())
    {
        auto& prop_vec = *it->second;
        auto const num_comp_mesh = prop_vec.getNumberOfGlobalComponents();
        if (num_comp_mesh == static_cast<int>(num_comp))
        {
            return prop_vec;
        }

        OGS_FATAL(
            "The requested property '{}' has {} components, but the one "
            "present in the mesh has {} components.",
            name, num_comp, num_comp_mesh);
    }

    auto const name_in_mesh = name + "_avg";
    auto [it, emplaced] = cell_averages_.emplace(
        name, MeshLib::getOrCreateMeshProperty<double>(
                  const_cast<MeshLib::Mesh&>(mesh_), name_in_mesh,
                  MeshLib::MeshItemType::Cell, num_comp));

    if (!it->second)
    {
        OGS_FATAL("The cell property '{}' could not be added to the mesh.",
                  name_in_mesh);
    }

    if (!emplaced)
    {
        OGS_FATAL(
            "Internal logic error. Something very bad happened. The cell "
            "property '{}' was not added to the list of cell averages to "
            "compute. There is some very strange inconsistency in the "
            "code. Trouble ahead!",
            name_in_mesh);
    }

    return *it->second;
}
}  // namespace ProcessLib
