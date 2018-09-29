/*!
  \file NodePartitionedMeshTester.cpp
  \author Wenqing Wang
  \date   2014.11
  \brief  Test class readNodePartitionedMesh to read node-wise partitioned mesh
  with MPI functions.

  \copyright
  Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <iomanip>
#include <fstream>
#include <string>

#include <mpi.h>

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/LogogCustomCout.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"

#include "MeshLib/IO/MPI_IO/NodePartitionedMeshReader.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

using namespace MeshLib;

int main(int argc, char *argv[])
{
    LOGOG_INITIALIZE();

    MPI_Init(&argc, &argv);

#ifdef USE_PETSC
    char help[] = "ogs6 with PETSc \n";
    PetscInitialize(&argc, &argv, nullptr, help);
#endif

    BaseLib::LogogCustomCout* out = new BaseLib::LogogCustomCout(1);
    using LogogFormatter = BaseLib::TemplateLogogFormatterSuppressedGCC
        <TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG>;
    LogogFormatter* fmt = new LogogFormatter();

    out->SetFormatter(*fmt);

    const std::string file_name = argv[1];
    std::string output_dir = "";
    if (argc > 2)
      output_dir = argv[2];

    NodePartitionedMesh *mesh = nullptr;
    {
        MeshLib::IO::NodePartitionedMeshReader read_pmesh(MPI_COMM_WORLD);
        mesh = read_pmesh.read(file_name);
    }
    if (!mesh)
    {
        ERR("Could not read mesh from files with prefix %s.", file_name.c_str());
        return EXIT_FAILURE;
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const std::string rank_str = std::to_string(rank);
    const std::string ofile_name = file_name + "_partition_" + rank_str + ".msh";
    std::ofstream os(BaseLib::joinPaths(output_dir, ofile_name), std::ios::trunc);

    // Output nodes
    os.setf(std::ios::scientific, std::ios::floatfield);
    std::setprecision(10);
    const std::size_t nn = mesh->getNumberOfNodes();
    for(std::size_t i=0; i<nn; i++)
    {
        const double *x = mesh->getNode(i)->getCoords();
        os << mesh->getNode(i)->getID() << " "
            << std::setw(14) << x[0]  << " " << x[1] << " "<< x[2] << "\n";
    }
    os.flush();

    // Output elements
    const std::size_t ne = mesh->getNumberOfElements();
    for(std::size_t i=0; i<ne; i++)
    {
        const Element *elem = mesh->getElement(i);
        Node* const* ele_nodes = elem->getNodes();

        for(unsigned j=0; j<elem->getNumberOfNodes(); j++)
        {
            os << ele_nodes[j]->getID() << " ";
        }
        os << "\n";
    }
    os.flush();

    delete mesh;

    delete out;
    delete fmt;

#ifdef USE_PETSC
    PetscFinalize();
#endif

    MPI_Finalize();

    LOGOG_SHUTDOWN();
}
