/*!
  \file NodePartitionedMeshReader.cpp
  \author Wenqing Wang
  \date   2014.08
  \brief  Define members of class NodePartitionedMeshReader to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodePartitionedMeshReader.h"

#include "logog/include/logog.hpp"

#include "BaseLib/RunTime.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"

using namespace MeshLib;

namespace FileIO
{
MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::read(MPI_Comm comm, const std::string &file_name_base)
{
    BaseLib::RunTime timer;
    timer.start();

    mpi_comm_ = comm;
    MPI_Comm_size(mpi_comm_, &_size);
    _size_str = std::to_string(_size);
    MPI_Comm_rank(mpi_comm_, &_rank);

    NodePartitionedMesh *mesh = nullptr;

    // Always try binary file first
    std::string fname_new = file_name_base + "_partitioned_msh_cfg" + _size_str + ".bin";

    if(!BaseLib::IsFileExisting(fname_new)) // doesn't exist binary file.
    {
        INFO("-->Reading ASCII mesh file ...");

        mesh = readASCII(file_name_base);
    }
    else
    {
        INFO("-->Reading binary mesh file ...");

        mesh = readBinary(file_name_base);
    }

    INFO("\t\n>>Total elapsed time in reading mesh:%f s\n", timer.elapsed());

    MPI_Barrier(comm);

    return mesh;
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader
::readBinary(const std::string &file_name_base)
{
    //----------------------------------------------------------------------------------
    // Read headers
    MPI_File fh;
    char ftype[] = "native";
    int file_status = 0;
    const std::string fname_header = file_name_base +  "_partitioned_msh_";
    const std::string fname_num_p_ext = _size_str + ".bin";

    std::string fname_new_base = fname_header + "cfg" + fname_num_p_ext;
    file_status = MPI_File_open(mpi_comm_, &fname_new_base[0], MPI_MODE_RDONLY,
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
    file_status = MPI_File_open(mpi_comm_, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status) // Failed to open the file
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    MPI_Datatype MPI_node;
    std::vector<NodeData> s_nodes(_num_nodes_part);
    buildNodeStrucTypeMPI(&s_nodes[0], &MPI_node);

    offset_new =  mesh_controls[10];
    MPI_File_set_view(fh, offset_new, MPI_node, MPI_node, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, &s_nodes[0], _num_nodes_part, MPI_node, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<unsigned> glb_node_ids;
    setNodes(s_nodes, mesh_nodes, glb_node_ids);

    MPI_Type_free(&MPI_node);

    //----------------------------------------------------------------------------------
    // Read non-ghost elements
    fname_new_base = fname_header +"ele" + fname_num_p_ext;
    file_status = MPI_File_open(mpi_comm_, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status) // Failed to open the file
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    long size_elem_info = _num_regular_elems_part + mesh_controls[8];
    std::vector<long> elem_data(size_elem_info);
    offset_new =  mesh_controls[11];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, &elem_data[0], size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    std::vector<MeshLib::Element*> mesh_elems(_num_regular_elems_part + _num_ghost_elems_part);
    setElements(mesh_nodes, elem_data.data(), mesh_elems);

    //----------------------------------------------------------------------------------
    //Read ghost element
    fname_new_base = fname_header + "ele_g" + fname_num_p_ext;
    file_status = MPI_File_open(mpi_comm_, &fname_new_base[0], MPI_MODE_RDONLY,
                                MPI_INFO_NULL, &fh);
    if(file_status)
    {
        printMessage(fname_new_base);
        return nullptr;
    }

    size_elem_info = _num_ghost_elems_part + mesh_controls[9];
    elem_data.resize(size_elem_info);
    offset_new =  mesh_controls[12];
    MPI_File_set_view(fh, offset_new, MPI_LONG, MPI_LONG, ftype, MPI_INFO_NULL);
    MPI_File_read(fh, &elem_data[0], size_elem_info, MPI_LONG, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    const bool process_ghost = true;
    setElements(mesh_nodes, elem_data.data(), mesh_elems, process_ghost);

    //----------------------------------------------------------------------------------
    //Create a mesh and return it
    const std::size_t n_nghost_elem = static_cast<std::size_t>(mesh_controls[2]);
    const unsigned n_global_base_nodes = static_cast<unsigned>(mesh_controls[6]);
    const unsigned n_global_nodes = static_cast<unsigned>(mesh_controls[7]);
    const unsigned n_base_nodes = static_cast<unsigned>(mesh_controls[1]);
    const unsigned n_active_base_nodes = static_cast<unsigned>(mesh_controls[4]);
    const unsigned n_active_nodes = static_cast<unsigned>(mesh_controls[5]);

    return new NodePartitionedMesh(BaseLib::extractBaseName(file_name_base) + _size_str,
                                   mesh_nodes, glb_node_ids,
                                   mesh_elems, n_nghost_elem,
                                   n_global_base_nodes, n_global_nodes,
                                   n_base_nodes, n_active_base_nodes, n_active_nodes);
}

bool NodePartitionedMeshReader::
openASCIIFiles(std::string const& file_name_base,
               std::ifstream& is_cfg, std::ifstream& is_node, std::ifstream& is_elem)
{
    const std::string fname_header = file_name_base +  "_partitioned_";
    const std::string fname_num_p_ext = _size_str + ".msh";

    std::string str_var = fname_header + "cfg"+ fname_num_p_ext;
    is_cfg.open(str_var.c_str());

    if( !is_cfg.good() )
    {
        printMessage(str_var);
        return false;
    }

    getline(is_cfg, str_var);
    int num_parts;
    is_cfg >> num_parts >> std::ws;

    if(num_parts != _size)
    {
        std::string str_var = "Aborting computation because of number of cores / subdomains mismatch.";
        printMessage(str_var, false);
        return false;
    }

    str_var = fname_header + "nodes_" + fname_num_p_ext;
    is_node.open(str_var.c_str());
    if( !is_node.good() )
    {
        return false;
    }

    str_var = fname_header + "elems_" + fname_num_p_ext;
    is_elem.open(str_var.c_str());
    if( !is_elem.good() )
    {
        printMessage(str_var);
        return false;
    }

    return true;
}

void NodePartitionedMeshReader
::readCastNodesASCII(std::ifstream& is_node, const int part_id,
                     std::vector<MeshLib::Node*> &mesh_nodes,
                     std::vector<unsigned> &glb_node_ids)
{
    int tag = 0;
    MPI_Status status;

    //----------------------------------------------------------------------------------
    // Read Nodes
    std::vector<NodeData> s_nodes(_num_nodes_part);

    MPI_Datatype MPI_node;
    if(part_id > 0)
        buildNodeStrucTypeMPI(&s_nodes[0], &MPI_node);

    if(_rank == 0)
    {
        for(long k=0; k<_num_nodes_part; k++)
        {
            NodeData *anode = &s_nodes[k];
            is_node >> anode->index
                    >> anode->x >> anode->y >> anode->z >> std::ws;
        }

        if(part_id == 0)
        {
            setNodes(s_nodes, mesh_nodes, glb_node_ids);
        }
        else
        {
            MPI_Send(&s_nodes[0], _num_nodes_part, MPI_node, part_id, tag, mpi_comm_);
        }
    }
    else if(part_id > 0 && _rank == part_id)
    {
        MPI_Recv(&s_nodes[0], _num_nodes_part, MPI_node, 0, tag, mpi_comm_, &status);
        setNodes(s_nodes, mesh_nodes, glb_node_ids);
    }

    if(part_id > 0)
        MPI_Type_free(&MPI_node);

}

void NodePartitionedMeshReader
::readCastElemsASCII(std::ifstream& is_elem, const int part_id,
                     const long data_size, const bool process_ghost,
                     const std::vector<MeshLib::Node*> &mesh_nodes,
                     std::vector<MeshLib::Element*> &mesh_elems)
{
    int tag = 0;
    MPI_Status status;

    long *elem_data = new long[data_size];
    if(_rank == 0)
    {
        readElementASCII(is_elem, elem_data, process_ghost);

        if(part_id == 0)
        {
            if(!process_ghost)
                mesh_elems.resize(_num_regular_elems_part + _num_ghost_elems_part);
            setElements(mesh_nodes, elem_data, mesh_elems, process_ghost);
        }
        else
        {
            MPI_Send(elem_data, data_size, MPI_LONG, part_id, tag, mpi_comm_);
        }
    }
    else if(part_id > 0 && _rank == part_id)
    {
        MPI_Recv(elem_data, data_size, MPI_LONG, 0, tag, mpi_comm_, &status);

        if(!process_ghost)
            mesh_elems.resize(_num_regular_elems_part + _num_ghost_elems_part);
        setElements(mesh_nodes, elem_data, mesh_elems, process_ghost);
    }

    delete [] elem_data;
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader
::readASCII(const std::string &file_name_base)
{
    std::ifstream is_cfg;
    std::ifstream is_node;
    std::ifstream is_elem;

    bool file_opened = false;
    if(_rank == 0)
    {
        file_opened = openASCIIFiles(file_name_base, is_cfg, is_node, is_elem);

        if(!file_opened)
        {
            is_cfg.close();
            is_node.close();
            is_elem.close();
        }
    }

    MPI_Bcast(&file_opened, 1, MPI_INT, 0, mpi_comm_);

    if(!file_opened)
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

    NodePartitionedMesh *np_mesh = nullptr;
    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<unsigned> glb_node_ids;
    std::vector<MeshLib::Element*> mesh_elems;

    for(int i=0; i<_size; i++)
    {
        if(_rank == 0)
        {
            INFO("-->Parallel reading the partitioned mesh: ");

            for(long j=0; j< num_controls; j++)
                is_cfg >> mesh_controls[j];
            is_cfg >> std::ws;
        }

        MPI_Bcast(mesh_controls, num_controls, MPI_LONG, 0, mpi_comm_);

        _num_nodes_part = mesh_controls[0];
        _num_regular_elems_part = mesh_controls[2];
        _num_ghost_elems_part = mesh_controls[3];

        //----------------------------------------------------------------------------------
        // Read Nodes
        readCastNodesASCII(is_node, i, mesh_nodes, glb_node_ids);

        //----------------------------------------------------------------------------------
        // Read elements
        bool process_ghost = false;
        const long size_elem_info = _num_regular_elems_part + mesh_controls[8];
        readCastElemsASCII(is_elem, i, size_elem_info, process_ghost,
                           mesh_nodes, mesh_elems);

        //-------------------------------------------------------------------------
        // Ghost elements
        process_ghost = true;
        const long size_elem_g_info = _num_ghost_elems_part + mesh_controls[9];
        readCastElemsASCII(is_elem, i, size_elem_g_info, process_ghost,
                           mesh_nodes, mesh_elems);

        if(_rank == i)
        {
            //----------------------------------------------------------------------------------
            //Create a mesh
            const std::size_t n_nghost_elem = static_cast<size_t>(mesh_controls[2]);
            const unsigned n_global_base_nodes = static_cast<unsigned>(mesh_controls[6]);
            const unsigned n_global_nodes = static_cast<unsigned>(mesh_controls[7]);
            const unsigned n_base_nodes = static_cast<unsigned>(mesh_controls[1]);
            const unsigned n_active_base_nodes = static_cast<unsigned>(mesh_controls[4]);
            const unsigned n_active_nodes = static_cast<unsigned>(mesh_controls[5]);

            np_mesh = new  NodePartitionedMesh(BaseLib::extractBaseName(file_name_base) + _size_str,
                                               mesh_nodes, glb_node_ids,
                                               mesh_elems, n_nghost_elem,
                                               n_global_base_nodes, n_global_nodes,
                                               n_base_nodes, n_active_base_nodes, n_active_nodes);
        }
    }

    if(_rank == 0)
    {
        is_cfg.close();
        is_node.close();
        is_elem.close();
    }

    MPI_Barrier(mpi_comm_);

    return  np_mesh;
}

void NodePartitionedMeshReader::readElementASCII(std::ifstream &ins,
        long *elem_data, const bool ghost)
{
    // Set number of elements.
    const long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    long counter = ne;
    for(long j=0; j<ne; j++)
    {
        elem_data[j] = counter;
        ins >> elem_data[counter];  //mat. idx
        counter++;
        ins >> elem_data[counter];  //type
        counter++;
        ins >> elem_data[counter];  //nnodes
        const long nn_e =  elem_data[counter];
        counter++;
        for(long k=0; k<nn_e; k++)
        {
            ins >> elem_data[counter];
            counter++;
        }
    }
}

void NodePartitionedMeshReader::setNodes(const std::vector<NodeData> &node_data,
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
        const long *elem_data, std::vector<MeshLib::Element*> &mesh_elems, const bool ghost)
{
    // Number of elements, ether ghost or regular
    const long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    const long id_offset = ghost ? _num_regular_elems_part : 0;
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

        mesh_elems[i + id_offset] = elem;
    }
}

void NodePartitionedMeshReader::buildNodeStrucTypeMPI(NodeData *anode, MPI_Datatype *mpi_node_ptr)
{
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

    MPI_Datatype my_comp_type[4];
    my_comp_type[0] = MPI_LONG;
    my_comp_type[1] = MPI_DOUBLE;
    my_comp_type[2] = MPI_DOUBLE;
    my_comp_type[3] = MPI_DOUBLE;

    // build datatype describing structure
    MPI_Type_create_struct(4, nblocklen, disp, my_comp_type, mpi_node_ptr);
    MPI_Type_commit(mpi_node_ptr);
}

void NodePartitionedMeshReader::printMessage(const std::string & err_message, const bool for_fileopen)
{
    if( for_fileopen )
    {
        ERR("! File %s does not exist.", err_message.c_str());
    }
    else
    {
        ERR( err_message.c_str() );
    }
}

}   // namespace FileIO
