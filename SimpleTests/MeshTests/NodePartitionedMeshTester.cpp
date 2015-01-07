/*!
  \file NodePartitionedMeshTester.cpp
  \author Wenqing Wang
  \date   2014.11
  \brief  Test class readNodePartitionedMesh to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#include <iomanip>
#include <fstream>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef USE_PETSC
#include <petscksp.h>
#endif

#include "logog/include/logog.hpp"

#include "BaseLib/LogogCustomCout.h"
#include "BaseLib/TemplateLogogFormatterSuppressedGCC.h"

#include "FileIO/MPI_IO/NodePartitionedMeshReader.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

using namespace MeshLib;
using namespace FileIO;

int main(int argc, char *argv[])
{
    LOGOG_INITIALIZE();
    {
#ifdef USE_MPI
        MPI_Init(&argc, &argv);
#endif

#ifdef USE_PETSC
        char help[] = "ogs6 with PETSc \n";
        PetscInitialize(&argc, &argv, nullptr, help);
#endif

        BaseLib::LogogCustomCout out(1);
        BaseLib::TemplateLogogFormatterSuppressedGCC<TOPIC_LEVEL_FLAG | TOPIC_FILE_NAME_FLAG | TOPIC_LINE_NUMBER_FLAG> custom_format;
        out.SetFormatter(custom_format);

        const std::string file_name = argv[1];

        NodePartitionedMeshReader read_pmesh;
        NodePartitionedMesh *mesh = read_pmesh.read(MPI_COMM_WORLD, file_name);

        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const std::string rank_str = std::to_string(rank);
        const std::string ofile_name = file_name + "_partition_" + rank_str + ".msh";
        std::ofstream os(ofile_name.data(), std::ios::trunc);

        // Output nodes
        os.setf(std::ios::scientific, std::ios::floatfield);
        std::setprecision(10);
        const size_t nn = mesh->getNNodes();
        for(size_t i=0; i<nn; i++)
        {
            const double *x = mesh->getNode(i)->getCoords();
            os << mesh->getNode(i)->getID() << " "
               << std::setw(14) << x[0]  << " " << x[1] << " "<< x[2] << "\n";
        }
        os.flush();

        // Output elements
        const size_t ne = mesh->getNElements();
        for(size_t i=0; i<ne; i++)
        {
            const Element *elem = mesh->getElement(i);
            Node* const* ele_nodes = elem->getNodes();

            for(unsigned j=0; j<elem->getNNodes(); j++)
            {
                os << ele_nodes[j]->getID() << " ";
            }
            os << "\n";
        }
        os.flush();

#ifdef USE_PETSC
        PetscFinalize();
#endif

#ifdef USE_MPI
        MPI_Finalize();
#endif

    } // make sure no logog objects exist when LOGOG_SHUTDOWN() is called.
    LOGOG_SHUTDOWN();
}

