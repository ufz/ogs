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
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

OGSMesh::OGSMesh(MeshLib::Mesh* mesh) : _mesh(mesh) {}

std::vector<double> OGSMesh::getPointCoordinates() const
{
    auto const& nodes = _mesh->getNodes();
    std::vector<double> coordinates;
    for (auto const& coords : nodes | MeshLib::views::coords)
    {
        std::copy(
            coords.begin(), coords.end(), std::back_inserter(coordinates));
    }
    return coordinates;
}
