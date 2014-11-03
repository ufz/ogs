/*!
  \file NodePartitionedMeshReader.cpp
  \author Wenqing Wang
  \date   2014.08
  \brief  Define members of class NodePartitionedMeshReader to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodePartitionedMeshReader.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "BaseLib/RunTime.h"
#include "BaseLib/FileTools.h"

#include "NodePartitionedMesh.h"
#include "Node.h"

#include "Elements/Element.h"
#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"

using namespace MeshLib;
using namespace std;

namespace FileIO
{
// Local function:
// Define MPI data type, MPI_Node_ptr, for struct MeshNode for palllel reading of nodes
static void buildNodeStrucTypeMPI(NodeData *anode, MPI_Datatype *MPI_Node_ptr)
{
    MPI_Datatype my_comp_type[4];
    my_comp_type[0] = MPI_LONG;
    my_comp_type[1] = MPI_DOUBLE;
    my_comp_type[2] = MPI_DOUBLE;
    my_comp_type[3] = MPI_DOUBLE;

    int nblocklen[4] = {1, 1, 1, 1};

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

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::read(MPI_Comm comm, const std::string &file_name)
{
    BaseLib::RunTime timer;
    timer.start();

    MPI_Comm_size(comm, &_size);
    _size_str = std::to_string(_size);
    MPI_Comm_rank(comm, &_rank);

    // Always try binary file first
    string fname_new = file_name + "_partitioned_msh_cfg" + _size_str + ".bin";

    NodePartitionedMesh *mesh = nullptr;

    MPI_File fh;
    int file_status = MPI_File_open(comm, &fname_new[0],
                                    MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    if(file_status)
    {
        if(_rank == 0)
            INFO("-->Reading ASCII mesh file ...");

        mesh = readASCII(comm, file_name);
    }
    else
    {
        MPI_File_close(&fh);
        if(_rank == 0)
            INFO("-->Reading binary mesh file ...");

        mesh = readBinary(comm, file_name);
    }

    if(_rank == 0)
        INFO("\t\n>>Total elapsed time in reading mesh:%f s\n", timer.elapsed());

    MPI_Barrier(comm);

    return mesh;
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::readBinary(MPI_Comm comm, const std::string &file_name)
{
    //----------------------------------------------------------------------------------
    // Read headers
    MPI_File fh;
    char ftype[] = "native";
    int file_status = 0;
    const string fname_header = file_name +  "_partitioned_msh_";
    const string fname_num_p_ext = _size_str + ".bin";

    string fname_new_base = fname_header + "cfg" + fname_num_p_ext;
    file_status = MPI_File_open(comm, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status) // Failed to open the file
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    const long num_controls = 14;
    /*
    Array that contains integer numbers of the partition data header
    and its size is 14.
            0:    Number of all nodes of a partition,
            1:    Number of nodes for linear element of a parition,
            2:    Number of non-ghost elements of a partition,
            3:    Number of ghost element of a partition,
            4:    Number of active nodes for linear element of a parition,
            5:    Number of all active nodes a parition,
            6:    Number of nodes for linear element of global mesh,
            7:    Number of all nodes of global mesh,
            8~12: Offsets of positions of partitions in the data arrays,
            13:   Reserved for exra flag.
    */
    long mesh_controls[num_controls];

    MPI_Offset offset_new;
    offset_new = _rank * num_controls * sizeof(long);
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, mesh_controls, num_controls, MPI_LONG, MPI_STATUS_IGNORE); //_all
    MPI_File_close(&fh);
    _num_nodes_part = mesh_controls[0];
    _num_regular_elems_part = mesh_controls[2];
    _num_ghost_elems_part = mesh_controls[3];

    //----------------------------------------------------------------------------------
    // Read Nodes
    fname_new_base = fname_header + "nod" + fname_num_p_ext;
    file_status = MPI_File_open(comm, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status) // Failed to open the file
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    MPI_Datatype MPI_node;
    NodeData *s_nodes = (NodeData *)malloc(sizeof(NodeData) * _num_nodes_part);
    buildNodeStrucTypeMPI(s_nodes, &MPI_node);

    offset_new =  mesh_controls[10];
    MPI_File_set_view(fh, offset_new, MPI_node, MPI_node, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, s_nodes, mesh_controls[0], MPI_node, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    vector<MeshLib::Node*> mesh_nodes;
    vector<unsigned> glb_node_ids;
    setNodes(s_nodes, mesh_nodes, glb_node_ids);

    free(s_nodes);
    MPI_Type_free(&MPI_node);

    //----------------------------------------------------------------------------------
    // Read non-ghost elements
    fname_new_base = fname_header +"ele" + fname_num_p_ext;
    file_status = MPI_File_open(comm, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status) // Failed to open the file
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    long size_elem_info = _num_regular_elems_part + mesh_controls[8];
    long *elem_data = (long *)malloc(sizeof(long) * size_elem_info );
    offset_new =  mesh_controls[11];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, elem_data, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::vector<MeshLib::Element*> mesh_elems;
    std::vector<short*> ghost_elems;
    setElements(mesh_nodes, elem_data, mesh_elems, ghost_elems);

    //----------------------------------------------------------------------------------
    //Read ghost element
    fname_new_base = fname_header + "ele_g" + fname_num_p_ext;
    file_status = MPI_File_open(comm, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status)
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    size_elem_info = _num_ghost_elems_part + mesh_controls[9];
    elem_data = (long *)realloc(elem_data, sizeof(long) * size_elem_info );
    offset_new =  mesh_controls[12];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, elem_data, size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    const bool process_ghost = true;
    setElements(mesh_nodes, elem_data, mesh_elems, ghost_elems, process_ghost);

    free(elem_data);

    //----------------------------------------------------------------------------------
    //Create a mesh and return it
    unsigned nnodes_active[] = { static_cast<unsigned>(mesh_controls[4]),
                                 static_cast<unsigned>(mesh_controls[5])
                               };
    unsigned nnodes_global[] = {static_cast<unsigned>(mesh_controls[6]),
                                static_cast<unsigned>(mesh_controls[7])
                               };

    return  new NodePartitionedMesh(BaseLib::extractBaseName(file_name) + _size_str, mesh_nodes,
                                    glb_node_ids, mesh_elems, ghost_elems,
                                    mesh_controls[2], nnodes_global, nnodes_active);
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::readASCII(MPI_Comm comm, const std::string &file_name)
{
    ifstream is_cfg;
    ifstream is_node;
    ifstream is_elem;
    int file_opened = 1;

    const string fname_header = file_name +  "_partitioned_";
    const string fname_num_p_ext = _size_str + ".msh";

    if(_rank == 0)
    {
        string str_var = fname_header + "cfg"+ fname_num_p_ext;
        is_cfg.open(str_var.c_str());

        if( !is_cfg.good() )
        {
            file_opened = false;
            printMessage(str_var);
        }
        else
        {
            getline(is_cfg, str_var);
            int num_parts;
            is_cfg >> num_parts >>ws;

            if(num_parts != _size)
            {
                file_opened = false;
                string str_var = "Aborting computation because of number of cores / subdomains mismatch.";
                printMessage(str_var, file_opened);
            }
        }

        str_var = fname_header + "nodes_" + fname_num_p_ext;
        is_node.open(str_var.c_str());
        if( !is_node.good() )
        {
            file_opened = false;
        }

        str_var = fname_header + "elems_" + fname_num_p_ext;
        is_elem.open(str_var.c_str());
        if( !is_elem.good() )
        {
            printMessage(str_var);
            file_opened = false;
        }
    }

    MPI_Bcast(&file_opened, 1, MPI_INT, 0, comm);

    if( !file_opened )
        return nullptr;

    const long num_controls = 11;
    /*
    Array that contains integer numbers of the partition data header
    and its size is 11.
            0:    Number of all nodes of a partition,
            1:    Number of nodes for linear element of a parition,
            2:    Number of non-ghost elements of a partition,
            3:    Number of ghost element of a partition,
            4:    Number of active nodes for linear element of a parition,
            5:    Number of all active nodes a parition,
            6:    Number of nodes for linear element of global mesh,
            7:    Number of all nodes of global mesh,
            8~9:  Offsets of positions of partitions in the data arrays,
            11:   Reserved for exra flag.
    */
    long mesh_controls[num_controls];

    NodePartitionedMesh * elem = nullptr;
    NodeData *s_nodes = nullptr;
    long *elem_info = nullptr;
    MPI_Datatype MPI_node;
    int tag[] = {0, 1, 2};
    //MPI_Request send_request, recv_request;
    MPI_Status status;
    vector<MeshLib::Node*> mesh_nodes;
    vector<unsigned> glb_node_ids;
    std::vector<MeshLib::Element*> mesh_elems;
    std::vector<short*> ghost_elems;
    for(int i=0; i<_size; i++)
    {
        if(_rank == 0)
        {
            cout<<"-->Parallel reading the partitioned mesh: "<<i<<endl;

            for(long j=0; j< num_controls; j++)
                is_cfg >> mesh_controls[j];
            is_cfg >> ws;

            if(i > 0)
            {
                MPI_Send(mesh_controls, num_controls, MPI_LONG, i, tag[0], comm);
            }
        }
        else if(i > 0)
        {
            if(_rank == i)
            {
                MPI_Recv(mesh_controls, num_controls, MPI_LONG, 0, tag[0], comm, &status);
            }
        }

        _num_nodes_part = mesh_controls[0];
        _num_regular_elems_part = mesh_controls[2];
        _num_ghost_elems_part = mesh_controls[3];

        //----------------------------------------------------------------------------------
        // Read Nodes
        s_nodes = (NodeData *)realloc(s_nodes, sizeof(NodeData) * mesh_controls[0]);

        if(i > 0)
            buildNodeStrucTypeMPI(s_nodes, &MPI_node);

        if(_rank == 0)
        {
            for(long k=0; k<_num_nodes_part; k++)
            {
                NodeData *anode = &s_nodes[k];
                is_node >> anode->index
                        >> anode->x >> anode->y >> anode->z >> ws;
            }

            if(i == 0)
            {
                setNodes(s_nodes, mesh_nodes, glb_node_ids);
            }
            else
            {
                MPI_Send(s_nodes, _num_nodes_part, MPI_node, i, tag[0], comm);
            }
        }
        else if(i > 0)
        {
            if(_rank == i)
            {
                MPI_Recv(s_nodes, _num_nodes_part, MPI_node, 0, tag[0], comm, &status);
                setNodes(s_nodes, mesh_nodes, glb_node_ids);
            }
        }

        if(i > 0)
            MPI_Type_free(&MPI_node);

        //----------------------------------------------------------------------------------
        // Read elements
        const long size_elem_info = _num_regular_elems_part + mesh_controls[8];
        elem_info = (long *)realloc(elem_info, sizeof(long) * size_elem_info );
        if(_rank == 0)
        {
            readElementASCII(is_elem, elem_info);
            if(i == 0)
            {
                setElements(mesh_nodes, elem_info, mesh_elems, ghost_elems);
            }
            else
            {
                MPI_Send(elem_info, size_elem_info, MPI_LONG, i, tag[1], comm);
            }
        }
        else if(i > 0)
        {
            if(_rank == i)
            {
                MPI_Recv(elem_info, size_elem_info, MPI_LONG, 0, tag[1], comm, &status);
                setElements(mesh_nodes, elem_info, mesh_elems, ghost_elems);
            }
        }

        //-------------------------------------------------------------------------
        // Ghost elements
        const bool process_ghost = true;
        const long size_elem_g_info = _num_ghost_elems_part + mesh_controls[9];
        elem_info = (long *)realloc(elem_info, sizeof(long) * size_elem_g_info );
        if(_rank == 0)
        {
            readElementASCII(is_elem, elem_info, true);

            if(i == 0)
            {
                setElements(mesh_nodes, elem_info, mesh_elems, ghost_elems, process_ghost);
            }
            else
            {
                MPI_Send(elem_info, size_elem_g_info, MPI_LONG, i, tag[2], comm);
            }
        }
        else if(i > 0)
        {
            if(_rank == i)
            {
                MPI_Recv(elem_info, size_elem_g_info, MPI_LONG, 0, tag[2], comm, &status);
                setElements(mesh_nodes, elem_info, mesh_elems, ghost_elems, process_ghost);
            }
        }

        if(_rank == i)
        {
            //----------------------------------------------------------------------------------
            //Create a mesh and return it
            unsigned nnodes_active[] = { static_cast<unsigned>(mesh_controls[4]),
                                         static_cast<unsigned>(mesh_controls[5])
                                       };
            unsigned nnodes_global[] = {static_cast<unsigned>(mesh_controls[6]),
                                        static_cast<unsigned>(mesh_controls[7])
                                       };

            elem =  new NodePartitionedMesh(BaseLib::extractBaseName(file_name) + _size_str, mesh_nodes,
                                            glb_node_ids, mesh_elems, ghost_elems,
                                            mesh_controls[2], nnodes_global, nnodes_active);
        }
    }

    if(s_nodes)
    {
        free(s_nodes);
    }
    if(elem_info)
    {
        free(elem_info);
    }

    if(_rank == 0)
    {
        is_cfg.close();
        is_node.close();
        is_elem.close();
    }

    MPI_Barrier(comm);

    return  elem;
}

void NodePartitionedMeshReader::readElementASCII(std::ifstream &ins,
        long *elem_info, const bool ghost)
{
    long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    long counter = ne;
    for(long j=0; j<ne; j++)
    {
        elem_info[j] = counter;
        ins >> elem_info[counter];  //mat. idx
        counter++;
        ins >> elem_info[counter];  //type
        counter++;
        ins >> elem_info[counter];  //nnodes
        const long nn_e =  elem_info[counter];
        counter++;
        for(long k=0; k<nn_e; k++)
        {
            ins >> elem_info[counter];
            counter++;
        }

        if( !ghost )
            continue;

        // ghost nodes for linear element
        ins >> elem_info[counter];
        counter++;

        ins >> elem_info[counter];
        const long nn_e_g =  elem_info[counter];
        counter++;
        for(long k=0; k<nn_e_g; k++)
        {
            ins >> elem_info[counter];
            counter++;
        }
    }
}

void NodePartitionedMeshReader::setNodes(const NodeData *node_data,
        std::vector<MeshLib::Node*> &mesh_node, std::vector<unsigned> &glb_node_ids)
{
    mesh_node.resize( _num_nodes_part );
    glb_node_ids.resize( _num_nodes_part );

    for(size_t i=0; i<mesh_node.size(); i++)
    {
        const NodeData *nd = &node_data[i];
        glb_node_ids[i] = static_cast<unsigned>(nd->index);
        mesh_node[i] = new MeshLib::Node(nd->x, nd->y, nd->z, i);
    }
}

void NodePartitionedMeshReader::setElements(const std::vector<MeshLib::Node*> &mesh_nodes,
        const long *elem_data, std::vector<MeshLib::Element*> &mesh_elems,
        std::vector<short*> &mesh_ghost_elems, const bool ghost)
{
    if( !ghost )
    {
        mesh_elems.resize(_num_regular_elems_part + _num_ghost_elems_part);
    }
    else
    {
        mesh_ghost_elems.resize(_num_ghost_elems_part);
    }

    // Number of elements, either ghost ot regular
    long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    long counter;
    for(long i=0; i<ne; i++)
    {
        counter = elem_data[i];

        const long mat_idx = static_cast<size_t>( elem_data[counter] );
        counter++;
        const long e_type = elem_data[counter];
        counter++;
        const long nnodes = elem_data[counter];
        counter++;

        MeshLib::Node **elem_nodes = new MeshLib::Node*[nnodes];
        for(long k=0; k<nnodes; k++)
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
            mesh_elems[i + _num_regular_elems_part] = elem;

            const short nn_gl = static_cast<short>( elem_data[counter] );
            counter++;
            const short nn_g = static_cast<short>( elem_data[counter] );
            counter++;

            short *local_ids = new short[nn_g+2];
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

void NodePartitionedMeshReader::printMessage(const std::string & err_message, const bool for_fileopen)
{
    if( for_fileopen )
    {
        if(_rank == 0)
            INFO("! File %s does not exist.", &err_message[0]);
    }
    else
    {
        if(_rank == 0)
            INFO( err_message.c_str() );
    }
}

} // end of name space FileIO
