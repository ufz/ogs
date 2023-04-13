/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

// Needs to be exported, see
// https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
class PYBIND11_EXPORT OGSMesh
{
public:
    explicit OGSMesh(MeshLib::Mesh* mesh);

    std::vector<double> getPointCoordinates() const;
    std::vector<double> getPointDataArray(std::string const& name) const;
    void setCellDataArray(std::string const& name,
                          std::vector<double> const& values);

private:
    MeshLib::Mesh* _mesh = nullptr;
};
