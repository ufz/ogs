/**
 * \file
 * \brief  Implementation of OpenGeoSys mesh python module
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
#include "OGSMesh.h"

/// python module name is OpenGeoSys
PYBIND11_MODULE(mesh, m)
{
    m.attr("__name__") = "ogs.mesh";
    m.doc() = "pybind11 ogs mesh example plugin";
    pybind11::class_<OGSMesh>(m, "OGSMesh")
        .def("getPointCoordinates", &OGSMesh::getPointCoordinates,
             "get node coordinates")
        .def("getCells", &OGSMesh::getCells,
             pybind11::return_value_policy::copy, "get cells")
        .def("getPointDataArray", &OGSMesh::getPointDataArray, "get point data")
        .def("setPointDataArray", &OGSMesh::setPointDataArray, "set point data")
        .def("setCellDataArray", &OGSMesh::setCellDataArray, "set cell data")
        .def("getCellDataArray", &OGSMesh::getCellDataArray, "get cell data");
}
