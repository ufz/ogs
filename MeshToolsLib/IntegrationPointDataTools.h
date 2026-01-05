// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <string>
#include <vector>

namespace MeshLib
{
class Element;
class Properties;
class PropertyVectorBase;
struct IntegrationPointMetaDataSingleField;
}  // namespace MeshLib

namespace MeshToolsLib
{
int getNumberOfElementIntegrationPoints(
    MeshLib::IntegrationPointMetaDataSingleField const& ip_meta_data,
    MeshLib::Element const& e);

std::vector<std::size_t> getIntegrationPointDataOffsetsOfMeshElements(
    std::vector<MeshLib::Element*> const& mesh_elements,
    MeshLib::PropertyVectorBase const& pv,
    MeshLib::Properties const& properties);
}  // namespace MeshToolsLib
