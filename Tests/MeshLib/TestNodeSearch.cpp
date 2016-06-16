/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "gtest/gtest.h"

#include <memory>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshGeoToolsLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshGenerators/RasterToMesh.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"

#include "GeoLib/Raster.h"

TEST(NodeSearch, UnusedNodes)
{
    std::array<double, 12> pix = {{0,0.1,0.2,0.1,0,0,0.1,0,0,0,-0.1,0}};
    GeoLib::RasterHeader const header =
        {4,3,MathLib::Point3d(std::array<double,3>{{0,0,0}}),1,-9999};
    GeoLib::Raster const raster(header,pix.begin(), pix.end());
    MeshLib::Mesh* mesh (MeshLib::RasterToMesh::convert(raster, MeshLib::MeshElemType::TRIANGLE, MeshLib::UseIntensityAs::ELEVATION));
    MeshLib::NodeSearch ns(*mesh);
    ns.searchUnused();
    std::vector<std::size_t> u_nodes = ns.getSearchedNodeIDs();
    ASSERT_EQ(0, u_nodes.size());

    std::vector<MeshLib::Node*> nodes = MeshLib::copyNodeVector(mesh->getNodes());
    nodes.push_back(new MeshLib::Node(-1,-1,-1));
    std::vector<MeshLib::Element*> elems = MeshLib::copyElementVector(mesh->getElements(),nodes);
    MeshLib::Mesh mesh2("mesh2", nodes, elems);
    MeshLib::NodeSearch ns2(mesh2);
    ns2.searchUnused();
    u_nodes = ns2.getSearchedNodeIDs();
    ASSERT_EQ(1, u_nodes.size());
    ASSERT_EQ(nodes.back()->getID(), u_nodes[0]);

    delete mesh;
}



