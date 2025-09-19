/**
 * \file
 * \brief  Implementation of OpenGeoSys mesh python module
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <numeric>

#include "BaseLib/ExportSymbol.h"
#include "BaseLib/Logging.h"
#include "InfoLib/GitInfo.h"
#include "OGSMesh.h"

/// python module name is OpenGeoSys
PYBIND11_MODULE(OGSMesh, m)
{
#ifndef NDEBUG
    BaseLib::initOGSLogger("all");
#else   // NDEBUG
    BaseLib::initOGSLogger("info");
#endif  // NDEBUG
    m.attr("__name__") = "ogs.mesh";
    m.doc() = "pybind11 ogs mesh example plugin";

    pybind11::enum_<MeshLib::MeshItemType>(m, "MeshItemType")
        .value("Node", MeshLib::MeshItemType::Node)
        .value("Edge", MeshLib::MeshItemType::Edge)
        .value("Face", MeshLib::MeshItemType::Face)
        .value("Cell", MeshLib::MeshItemType::Cell)
        .value("IntegrationPoint", MeshLib::MeshItemType::IntegrationPoint)
        .export_values()
        // Add nice string conversion
        .def("__str__",
             [](MeshLib::MeshItemType t) { return std::string(toString(t)); });

    pybind11::class_<OGSMesh>(m, "OGSMesh")
        .def("points", &OGSMesh::getPointCoordinates, "get node coordinates")
        .def("cells", &OGSMesh::getCells, pybind11::return_value_policy::copy,
             "get cells")
        .def("data_array_names", &OGSMesh::getDataArrayNames,
             "get names of all data arrays / property "
             "vectors stored in the mesh")
        .def("mesh_item_type", &OGSMesh::meshItemType, pybind11::arg("name"),
             "returns MeshItemType")
        .def("data_array", &OGSMesh::dataArray_dispatch,
             pybind11::return_value_policy::reference, pybind11::arg("name"),
             pybind11::arg("dtype"),
             "Access a data array / property vector (accesses OGS memory "
             "directly using numpy array with appropriate shape)");
}
