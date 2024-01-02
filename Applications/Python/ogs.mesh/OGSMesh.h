/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <pybind11/pybind11.h>

#include <string>
#include <utility>
#include <vector>

namespace MeshLib
{
class Mesh;
}

// Needs to be exported, see
// https://pybind11.readthedocs.io/en/stable/advanced/misc.html#partitioning-code-over-multiple-extension-modules
class PYBIND11_EXPORT OGSMesh
{
public:
    explicit OGSMesh(MeshLib::Mesh& mesh);

    std::vector<double> getPointCoordinates() const;
    void setPointDataArray(std::string const& name,
                           std::vector<double> const& values,
                           std::size_t const number_of_components);
    std::vector<double> getPointDataArray(
        std::string const& name,
        std::size_t const number_of_components = 1) const;
    std::pair<std::vector<int>, std::vector<int>> getCells() const;
    void setCellDataArray(std::string const& name,
                          std::vector<double> const& values,
                          std::size_t const number_of_components);
    std::vector<double> getCellDataArray(
        std::string const& name, std::size_t const number_of_components) const;

private:
    MeshLib::Mesh& _mesh;
};
