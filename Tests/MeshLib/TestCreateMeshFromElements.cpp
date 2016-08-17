/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>
#include <cmath>

#include "gtest/gtest.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"

#include "MeshLib/Utils/createMeshFromElements.h"

// Generates new meshes from subsets of elements of a regular hex mesh
TEST(MeshLib, createMeshFromElements)
{
    unsigned n_x(10);
    unsigned n_y(5);
    unsigned n_z(2);
    double delta(1.0);
    std::unique_ptr<MeshLib::Mesh> hex_mesh(
        MeshLib::MeshGenerator::generateRegularHexMesh(n_x, n_y, n_z, delta));

    std::vector<MeshLib::Element*> subset_elements(hex_mesh->getElements());

    // use all elements of the hexahedron mesh to create a new mesh
    std::unique_ptr<MeshLib::Mesh> subset_mesh(
        MeshLib::createMeshFromElements(subset_elements, "CompleteMesh"));
    ASSERT_EQ(hex_mesh->getNumberOfElements(),
              subset_mesh->getNumberOfElements());
    ASSERT_EQ(hex_mesh->getNumberOfNodes(),
              subset_mesh->getNumberOfNodes());

    // remove one element of the hexahedron mesh and create a new mesh from the
    // remaining elements
    subset_elements.pop_back();
    subset_mesh.reset(
        MeshLib::createMeshFromElements(subset_elements, "OneElementRemoved"));
    ASSERT_EQ(hex_mesh->getNumberOfElements()-1,
              subset_mesh->getNumberOfElements());
    ASSERT_EQ(hex_mesh->getNumberOfNodes()-1,
              subset_mesh->getNumberOfNodes());

    // remove complete layer of the hexahedron mesh and create a new mesh from
    // the remaining elements
    for (std::size_t k(0); k<(n_x*n_y-1); ++k)
        subset_elements.pop_back();

    subset_mesh.reset(
        MeshLib::createMeshFromElements(subset_elements, "OneLayerRemoved"));
    ASSERT_EQ(n_x*n_y*(n_z-1),
              subset_mesh->getNumberOfElements());
    ASSERT_EQ((n_x+1)*(n_y+1)*n_z,
              subset_mesh->getNumberOfNodes());
}
