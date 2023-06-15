/**
 * \file
 * \brief  Implementation of OGSMesh
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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

#include <numeric>
#include <range/v3/numeric.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

#include "BaseLib/ExportSymbol.h"
#include "BaseLib/Logging.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "MeshLib/VtkOGSEnum.h"

OGSMesh::OGSMesh(MeshLib::Mesh& mesh) : _mesh(mesh) {}

std::vector<double> OGSMesh::getPointCoordinates() const
{
    auto const& nodes = _mesh.getNodes();
    std::vector<double> coordinates;
    for (auto const& coords : nodes | MeshLib::views::coords)
    {
        std::copy(
            coords.begin(), coords.end(), std::back_inserter(coordinates));
    }
    return coordinates;
}

std::pair<std::vector<int>, std::vector<int>> OGSMesh::getCells() const
{
    auto const& elements = _mesh.getElements();
    std::vector<int> cells;
    std::vector<int> cell_types;
    for (auto const* element : elements)
    {
        auto const number_of_nodes =
            static_cast<int>(element->getNumberOfNodes());
        cells.push_back(number_of_nodes);
        for (int i = 0; i < number_of_nodes; ++i)
        {
            cells.push_back(element->getNode(i)->getID());
        }
        cell_types.push_back(OGSToVtkCellType(element->getCellType()));
    }
    return {cells, cell_types};
}

void OGSMesh::setPointDataArray(std::string const& name,
                                std::vector<double> const& values,
                                std::size_t const number_of_components)
{
    auto* pv = MeshLib::getOrCreateMeshProperty<double>(
        _mesh, name, MeshLib::MeshItemType::Node, number_of_components);
    if (pv == nullptr)
    {
        OGS_FATAL("Couldn't access cell/element property '{}'.", name);
    }
    if (pv->size() != values.size())
    {
        OGS_FATAL(
            "OGSMesh::setPointDataArray: size mismatch: property vector has "
            "size '{}', while the number of values is '{}'.",
            pv->size(), values.size());
    }
    std::copy(values.begin(), values.end(), pv->data());
}

std::vector<double> OGSMesh::getPointDataArray(
    std::string const& name, std::size_t const number_of_components) const
{
    auto const* pv = MeshLib::getOrCreateMeshProperty<double>(
        _mesh, name, MeshLib::MeshItemType::Node, number_of_components);
    if (pv == nullptr)
    {
        OGS_FATAL("Couldn't access point/node property '{}'.", name);
    }
    std::vector<double> data_array;
    data_array.reserve(pv->getNumberOfTuples() * number_of_components);
    std::copy(pv->begin(), pv->end(), std::back_inserter(data_array));

    return data_array;
}

void OGSMesh::setCellDataArray(std::string const& name,
                               std::vector<double> const& values,
                               std::size_t const number_of_components)
{
    auto* pv = MeshLib::getOrCreateMeshProperty<double>(
        _mesh, name, MeshLib::MeshItemType::Cell, number_of_components);
    if (pv == nullptr)
    {
        OGS_FATAL("Couldn't access cell/element property '{}'.", name);
    }
    if (pv->size() != values.size())
    {
        OGS_FATAL(
            "OGSMesh::setCellDataArray: size mismatch: property vector has "
            "size '{}', while the number of values is '{}'.",
            pv->size(), values.size());
    }
    std::copy(values.begin(), values.end(), pv->data());
}

std::vector<double> OGSMesh::getCellDataArray(
    std::string const& name, std::size_t const number_of_components) const
{
    auto const* pv = MeshLib::getOrCreateMeshProperty<double>(
        _mesh, name, MeshLib::MeshItemType::Cell, number_of_components);
    if (pv == nullptr)
    {
        OGS_FATAL("Couldn't access cell/element property '{}'.", name);
    }
    std::vector<double> data_array;
    data_array.reserve(pv->getNumberOfTuples() * number_of_components);
    std::copy(pv->begin(), pv->end(), std::back_inserter(data_array));

    return data_array;
}
