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

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include <cstdlib>

#include "BaseLib/StringTools.h" // for int to string
#include "BaseLib/WallClockTimer.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/NodePartitionedMes.h"
#include "MeshLib/Node.h"

#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"

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
            INFO("-->Reading ASCII mesh file ...");


        mesh = readASCII(comm, file_name);
    }
    else
    {
        MPI_File_close(&fh);
        if( _rank == 0 )
            INFO("-->Reading binary mesh file ...");

        mesh = readBinary(comm, file_name);
    }

    if( _rank == 0 )
        INFO("\t\n>>Total elapsed time in reading mesh:%f s\n", timer.elapsed());


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
    if (rc ) // Failed to open the file
    {
        quit(fname_new_base);
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
    if (rc ) // Failed to open the file
    {
        quit(fname_new_base);
    }

    MPI_Datatype MPI_node;
    NodeData *s_nodes = (NodeData *)malloc(sizeof(NodeData) * _mesh_controls[0]);
    buildNodeStrucTypeMPI(s_nodes, &MPI_node);

    offset_new =  _mesh_controls[10];
    MPI_File_set_view(fh, offset_new, MPI_node, MPI_node,  &ftype[0], MPI_INFO_NULL);
    MPI_File_read(fh, s_nodes, _mesh_controls[0], MPI_node, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    vector<MeshLib::Node*> mesh_nodes;
    setNodes(s_nodes, mesh_nodes);
    free(s_nodes);

    MPI_Type_free(&MPI_node);

    // Elements
    fname_new_base = file_name + "_partitioned_msh_ele" + _size_str + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if (rc ) // Failed to open the file
    {
        quit(fname_new_base);
    }

    MyInt size_elem_info =   _mesh_controls[2] + _mesh_controls[8];
    MyInt *elem_data = (MyInt *)malloc(sizeof(MyInt) * size_elem_info );
    offset_new =  _mesh_controls[11];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, &ftype[0], MPI_INFO_NULL);
    MPI_File_read(fh, elem_data, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    // _mesh_controls[2]: number of regular elements. _mesh_controls[3]: number of ghost elements
    std::vector<MeshLib::Element*> mesh_elems(_mesh_controls[2] + _mesh_controls[3]);
    std::vector<short*> ghost_elems(0);
    setElements(mesh_nodes, elem_data, mesh_elems, ghost_elems);

    //Ghost element
    fname_new_base = file_name + "_partitioned_msh_ele_g" + _size_str + ".bin";
    rc = MPI_File_open(MPI_COMM_WORLD, &fname_new_base[0], MPI_MODE_RDONLY,
                       MPI_INFO_NULL, &fh);
    if( rc )
    {
        MPI_Finalize();
        if( _rank == 0 )
            INFO("! File %s does not exist.", &fname_new_base[0]);
        exit(EXIT_FAILURE);
    }

    size_elem_info =   _mesh_controls[3] + _mesh_controls[9];
    elem_data = (MyInt *)realloc(elem_data, sizeof(MyInt) * size_elem_info );
    offset_new =  _mesh_controls[12];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG,  &ftype[0], MPI_INFO_NULL);
    MPI_File_read(fh, elem_data, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::vector<short*> ghost_elems( _mesh_controls[3] );
    setElements(mesh_nodes, elem_data, mesh_elems, ghost_elems);

    free(elem_info);

    unsigned nnodes_local[] = { _mesh_controls[4], _mesh_controls[5] };
    unsigned nnodes_global[] = { _mesh_controls[6], _mesh_controls[7] };

    return  new NodePartitionedMesh(file_name + _size_str,
                                    mesh_nodes,  mesh_elems, nnodes_global, nnodes_local);
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
        mesh_node[i] = new MeshLib::Node( nd->x, nd->y, nd->z, static_cast<unsigned>(nd->index) );
    }
}

void readNodePartitionedMesh::setElements(const std::vector<MeshLib::Node*> &mesh_nodes,
        const MyInt *elem_data, std::vector<MeshLib::Element*> &mesh_elems,
        std::vector<unsigned*> &mesh_ghost_elems)
{
    bool ghost = (mesh_ghost_elems.size() > 0) ? true : false;
    // Number of elements, either ghost ot regular
    MyInt ne = ghost ? _mesh_controls[3] : _mesh_controls[2];
    MyInt counter;
    for(MyInt i=0; i<ne; i++)
    {
        counter = elem_data[i];

        const MyInt mat_idx = static_cast<size_t>( elem_data[counter] );
        counter++;
        const MyInt e_type = elem_data[counter];
        counter++;
        const MyInt nnodes = elem_data[counter];
        counter++;

        MeshLib::Node **elem_nodes = new MeshLib::Node*[nnodes]
        for(MyInt k=0; k<nnodes; k++)
        {
            elem_nodes[k] = mesh_nodes[ elem_data[counter] ];
            counter++;
        }

        MeshLib::Element *elem = nullptr;

        switch(e_type)
        {
            case 1:
                elem = new MeshLib::Line(elem_nodes, mat_idx);
                break;
            case 2:
                elem = new MeshLib::Quad(elem_nodes, mat_idx);
                break;
            case 3:
                elem = new MeshLib::Hex(elem_nodes, mat_idx);
                break;
            case 4:
                elem = new MeshLib::Tri(elem_nodes, mat_idx);
                break;
            case 5:
                elem = new MeshLib::Tet(elem_nodes, mat_idx);
                break;
            case 6:
                elem = new MeshLib::Prism(elem_nodes, mat_idx);
                break;
            case 7:
                elem = new MeshLib::Pyramid(elem_nodes, mat_idx);
                break;
        }

        if( ghost )
        {
            mesh_elems[i + _mesh_controls[2] ] = elem;

            const short nn_gl = static_cast<short>( elem_data[counter] );
            counter++;
            const short nn_g = static_cast<short>( elem_data[counter] );
            counter++;

            short *local_ids = new int[nn_g+2];
            local_ids[0] = nn_gl;
            local_ids[1] = nn_g;
            for(int k=2; k<nn_g+2; k++)
            {
                local_ids[k] = static_cast<short>( elem_data[counter] );
                counter++;
            }

            mesh_ghost_elems[i] = local_ids;
        }
        else
        {
            mesh_elems[i] = elem;
        }

    }
}

void readNodePartitionedMesh::quit(const string & file_name)
{
    MPI_Finalize();
    if(_rank == 0 )
        INFO("! File %s does not exist.", &file_name[0]);
    exit(EXIT_FAILURE);
}
} // end of name space FileIO
