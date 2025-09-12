/**
 * \file
 * \brief  Implementation of OGSMesh
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "OGSMesh.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <spdlog/spdlog.h>

#include <range/v3/view/join.hpp>
#include <vector>

#include "BaseLib/ExportSymbol.h"
#include "BaseLib/Logging.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/VtkOGSEnum.h"

OGSMesh::OGSMesh(MeshLib::Mesh& mesh) : _mesh(mesh) {}

std::vector<double> OGSMesh::getPointCoordinates() const
{
    return ranges::to<std::vector>(_mesh.getNodes() | MeshLib::views::coords |
                                   ranges::views::join);
}

std::pair<std::vector<int>, std::vector<int>> OGSMesh::getCells() const
{
    auto const& elements = _mesh.getElements();
    std::vector<int> cells;
    std::vector<int> cell_types;
    for (auto const* element : elements)
    {
        ranges::copy(element->nodes() | MeshLib::views::ids,
                     std::back_inserter(cells));
        cell_types.push_back(OGSToVtkCellType(element->getCellType()));
    }
    return {cells, cell_types};
}

std::vector<std::string> OGSMesh::getDataArrayNames() const
{
    return _mesh.getProperties().getPropertyVectorNames();
}

MeshLib::MeshItemType OGSMesh::meshItemType(std::string_view const name) const
{
    auto const& properties = _mesh.getProperties();

    auto const& found = BaseLib::findElementOrError(
        properties,
        [&name](auto const& p) { return p.first == name; },
        [&name]()
        {
            OGS_FATAL("A property with the name '{}' doesn't exist in the mesh",
                      name);
        });
    return found.second->getMeshItemType();
}
