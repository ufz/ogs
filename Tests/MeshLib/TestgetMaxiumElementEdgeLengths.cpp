/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 13, 2022, 10:40 AM
 */

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <memory>
#include <vector>

#include "InfoLib/TestInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getMaxiumElementEdgeLengths.h"

TEST(MeshLib, getMaxiumElementEdgeLengths)
{
    // Use mesh quadratic_mesh.vtu that contains 8 elements with different
    // types:
    std::string const file_name =
        TestInfoLib::TestInfo::data_path + "/Utils/GMSH2OGS/quadratic_mesh.vtu";
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(file_name));
    if (!mesh)
    {
        OGS_FATAL("Could not read mesh from '{:s}' file. No mesh created.",
                  file_name);
    }

    auto const element_sizes =
        MeshLib::getMaxiumElementEdgeLengths(mesh->getElements());
    std::vector<double> const expected_element_sizes = {
        1, std::sqrt(2), 1, 1, std::sqrt(2), 1, std::sqrt(2), std::sqrt(1.5)};

    for (std::size_t i = 0; i < element_sizes.size(); i++)
    {
        ASSERT_NEAR(expected_element_sizes[i],
                    element_sizes[i],
                    std::numeric_limits<double>::epsilon());
    }
}
