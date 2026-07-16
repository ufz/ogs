// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "ParameterLib/Parameter.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshToolsLib
{
struct InitialConditionDataSet
{
    std::string variable_name;
    std::vector<std::size_t> element_ids_for_selected_materials;
    std::unique_ptr<const MeshLib::PropertyVector<double>> property_copy;
    MeshLib::Mesh* mesh = nullptr;
    ParameterLib::Parameter<double> const* parameter = nullptr;
    MeshLib::MeshItemType mesh_item_type = MeshLib::MeshItemType::Cell;
};

void overwriteMeshFieldDataByMaterialIDs(
    std::vector<InitialConditionDataSet>& data_initial_conditions);
}  // namespace MeshToolsLib
