/*!
  \file TestNodePartitionedMeshReader.cpp
  \author Wenqing Wang
  \date   2014.08
  \brief  Test class readNodePartitionedMesh to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/


#include "gtest/gtest.h"
#include "BaseLib/BuildInfo.h" // Path to input file

#include "Node.h"

#include "MPI_MeshIO/NodePartitionedMeshReader.h"

using namespace MeshLib;
using namespace FileIO;

TEST(MPITest_FileIO, TestNodePartitionedMeshReader)
{
    std::string fname = BaseLib::BuildInfo::source_path + "/Tests/FileIO/PartitionedMesh/mesh_3d";

    NodePartitionedMeshReader read_pmesh;

    NodePartitionedMesh *mesh = read_pmesh.read(MPI_COMM_WORLD, fname);


    const size_t nn = mesh->getNNodes();
    const size_t ne = mesh->getNElements();
    const unsigned id_1stg_elem = mesh->getStartIndexOfGhostElement();
    const unsigned ne_ghost = ne - id_1stg_elem;

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch(rank)
    {
        case 0:
            {
                ASSERT_EQ(42u, nn);
                ASSERT_EQ(98u, ne);
                ASSERT_EQ(35u, id_1stg_elem);

                const double *x0 = mesh->getNode(0)->getCoords();
                ASSERT_NEAR(0.0, x0[0], 1.e-10);
                ASSERT_NEAR(0.0, x0[1], 1.e-10);
                ASSERT_NEAR(1.0, x0[2], 1.e-10);
                const double *x1 = mesh->getNode(nn-1)->getCoords();
                ASSERT_NEAR(4.99996114718444e-01, x1[0], 1.e-16);
                ASSERT_NEAR(7.50004016082884e-01, x1[1], 1.e-16);
                ASSERT_NEAR(2.50030875593138e-01, x1[2], 1.e-16);

                unsigned e_act_nn = mesh->getElementActiveNNodes(ne_ghost-1);
                ASSERT_EQ(1u, e_act_nn);
                short *e_act_nodes = mesh->getElementActiveNodes(ne_ghost-1);
                ASSERT_EQ(3u, e_act_nodes[0]);
            }
            break;
        case 1:
            {
                ASSERT_EQ(44u, nn);
                ASSERT_EQ(104u, ne);
                ASSERT_EQ(31u, id_1stg_elem);

                const double *x0 = mesh->getNode(0)->getCoords();
                ASSERT_NEAR(0.0, x0[0], 1.e-10);
                ASSERT_NEAR(0.0, x0[1], 1.e-10);
                ASSERT_NEAR(0.0, x0[2], 1.e-10);
                const double *x1 = mesh->getNode(nn-1)->getCoords();
                ASSERT_NEAR(2.50078531617154e-01, x1[0], 1.e-16);
                ASSERT_NEAR(2.50077599298029e-01, x1[1], 1.e-16);
                ASSERT_NEAR(1.0, x1[2], 1.e-16);

                unsigned e_act_nn = mesh->getElementActiveNNodes(ne_ghost-1);
                ASSERT_EQ(1u, e_act_nn);
                short *e_act_nodes = mesh->getElementActiveNodes(ne_ghost-1);
                ASSERT_EQ(1u, e_act_nodes[0]);
            }
            break;
        case 2:
            {
                ASSERT_EQ(42u, nn);
                ASSERT_EQ(96u, ne);
                ASSERT_EQ(30u, id_1stg_elem);

                const double *x0 = mesh->getNode(0)->getCoords();
                ASSERT_NEAR(1.0, x0[0], 1.e-10);
                ASSERT_NEAR(0.0, x0[1], 1.e-10);
                ASSERT_NEAR(0.0, x0[2], 1.e-10);
                const double *x1 = mesh->getNode(nn-1)->getCoords();
                ASSERT_NEAR(5.00082423503819e-01, x1[0], 1.e-16);
                ASSERT_NEAR(2.50043690294878e-01, x1[1], 1.e-16);
                ASSERT_NEAR(7.50023794732280e-01, x1[2], 1.e-16);

                unsigned e_act_nn = mesh->getElementActiveNNodes(ne_ghost-1);
                ASSERT_EQ(1u, e_act_nn);
                short *e_act_nodes = mesh->getElementActiveNodes(ne_ghost-1);
                ASSERT_EQ(0u, e_act_nodes[0]);
            }
            break;
        default:
            break;
    }

    delete mesh;
}
