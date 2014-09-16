/*!
  \file readNodePartitionedMesh.cpp
  \author Wenqing Wang
  \date   2014.08
  \brief  Define members of class readNodePartitionedMesh to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "readNodePartitionedMesh.h"

#include <cstdlib>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/NodePartitionedMes.h"
#include "MeshLib/Node.h"

#include "BaseLib/StringTools.h" // for int to string
#include "BaseLib/WallClockTimer.h"

using namespace MeshLib;

namespace FileIO
{
void buildNodeStrucTypeMPI(NodeData *anode, MPI_Datatype *MPI_Node_ptr)
{
    MPI_Datatype my_comp_type[4];
    my_comp_type[0] = MPI_LONG;
    my_comp_type[1] = MPI_DOUBLE;
    my_comp_type[2] = MPI_DOUBLE;
    my_comp_type[3] = MPI_DOUBLE;

    int nblocklen[4];
    nblocklen[0] = 1;
    nblocklen[1] = 1;
    nblocklen[2] = 1;
    nblocklen[3] = 1;

    MPI_Aint disp[4], base;
    MPI_Get_address(anode, disp);
    MPI_Get_address(&(anode[0].x), disp+1);
    MPI_Get_address(&(anode[0].y), disp+2);
    MPI_Get_address(&(anode[0].z), disp+3);
    base = disp[0];
    for(int j=0; j <4; j++)
    {
        disp[j] -= base;
    }

    // build datatype describing structure
    MPI_Type_create_struct(4, nblocklen, disp, my_comp_type, MPI_Node_ptr);
    MPI_Type_commit(MPI_Node_ptr);
}

MeshLib::NodePartitionedMesh* readNodePartitionedMesh::read(MPI_Comm comm, const std::string &file_name)
{
    WallClockTimer timer;
    timer.start();

    MPI_Comm_size(comm, &_size);
    string _size_str = number2str(_size);
    MPI_Comm_rank(comm, &_rank);

    // Always try binary file first
    string fname_new = file_name + "_partitioned_msh_cfg" + _size_str + ".bin";

    NodePartitionedMesh *mesh = nullptr;

    MPI_File fh;
    int rc = 0;
    rc = MPI_File_open(comm, &fname_new[0],
                       MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if (rc )
    {

        if( _rank == 0 )
            printf("-->Reading ASCII mesh file ...");

        mesh = readASCII(comm, file_name);
    }
    else
    {
        MPI_File_close(&fh);
        if( _rank == 0 )
            printf("-->Reading binary mesh file ...");

        mesh = readBinary(comm, file_name);
    }

    if( _rank == 0 )
        printf( "\t\n>>Total elapsed time in reading mesh:%f s\n", timer.elapsed() );

    MPI_Barrier(comm);

    return mesh;
}

MeshLib::NodePartitionedMesh* readNodePartitionedMesh::readBinary(MPI_Comm comm, const std::string &file_name);
{
    // Read headers
    MPI_File fh;
    string ftype = "native";
    int rc = 0;

    string fname_new_base = file_name + "_partitioned_msh_cfg" + _size_str + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc )
    {

        MPI_Finalize();
        if( _rank == 0 )
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit( EXIT_FAILURE );
    }
    //
    MPI_Offset offset_new;
    offset_new = _rank * _mesh_controls * sizeof(MyInt);
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG,  &ftype[0], MPI_INFO_NULL);
    MPI_File_read(fh, _mesh_controls, _num_controls, MPI_LONG, MPI_STATUS_IGNORE); //_all
    MPI_File_close(&fh);

    // Nodes
    fname_new_base = file_name + "_partitioned_msh_nod" + _size_str + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc )
    {

        MPI_Finalize();
        if( _rank == 0 )
            printf("! File %s does not exist.", &fname_new_base[0]);
        exit(EXIT_FAILURE);
    }

    MPI_Datatype MPI_node;
    NodeData *s_nodes = (NodeData *)malloc(s_nodes, sizeof(NodeData) * _mesh_controls[0]);
    buildNodeStrucTypeMPI(s_nodes, &MPI_node);

    offset_new =  _mesh_controls[10];
    MPI_File_set_view(fh, offset_new, MPI_node, MPI_node,  &ftype[0], MPI_INFO_NULL);
    MPI_File_read(fh, s_nodes, _mesh_controls[0], MPI_node, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    mesh->setSubdomainNodes(mesh_header, s_nodes);
    free(s_nodes);

    MPI_Type_free(&MPI_node);
}

MeshLib::NodePartitionedMesh* readNodePartitionedMesh::readASCII(MPI_Comm comm, const std::string &file_name)
{

}

void readNodePartitionedMesh::setNodes(const NodeData *node_data,
                                       std::vector<MeshLib::Node*> &mesh_node)
{
    mesh_node.resize( _mesh_controls[0] );

    for(size_t i=0; i<mesh_node.size(); i++)
    {
        const NodeData *nd = &node_data[i];
        mesh_node[i] = new MeshLib::Node(nd->x, nd->y, nd->z, nd->index);
    }
}

void readNodePartitionedMesh::setElements(const MyInt *elem_data,
        std::vector<MeshLib::Element*> &mesh_elem,
        const bool ghost)
{
}

} // end of name space FileIO
