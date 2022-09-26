/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on September 26, 2022, 11:32 AM
 */

#pragma once
#include <string>
#include <vector>

namespace MeshLib
{
class Element;
class Properties;
class PropertyVectorBase;
}  // namespace MeshLib

namespace ApplicationUtils
{

int getNumberOfElementIntegrationPoints(std::string const& data_name,
                                        MeshLib::Properties const& properties,
                                        MeshLib::Element const& e);

std::vector<std::size_t> getIntegrationPointDataOffsetsOfMeshElements(
    std::vector<MeshLib::Element*> const& mesh_elements,
    MeshLib::PropertyVectorBase const& pv,
    MeshLib::Properties const& properties);
}  // namespace ApplicationUtils
