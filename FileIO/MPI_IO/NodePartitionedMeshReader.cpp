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

#include <array>

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

namespace FileIO
{

NodePartitionedMeshReader::NodePartitionedMeshReader(MPI_Comm comm)
    : _mpi_comm(comm)
{
    MPI_Comm_size(_mpi_comm, &_mpi_comm_size);
    MPI_Comm_rank(_mpi_comm, &_mpi_rank);

    registerNodeDataMpiType();
}

NodePartitionedMeshReader::~NodePartitionedMeshReader()
{
    MPI_Type_free(&_mpi_node_type);
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::read(
    const std::string &file_name_base)
{
    BaseLib::RunTime timer;
    timer.start();

    MeshLib::NodePartitionedMesh *mesh = nullptr;

    // Always try binary file first
    std::string const fname_new = file_name_base + "_partitioned_msh_cfg" +
        std::to_string(_mpi_comm_size) + ".bin";

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

    MPI_Barrier(_mpi_comm);

    return mesh;
}

template <typename DATA>
bool
NodePartitionedMeshReader::readBinaryDataFromFile(std::string const& filename,
        MPI_Offset offset, MPI_Datatype type, DATA& data) const
{
    // Open file
    MPI_File file;

    char* filename_char = const_cast<char*>(filename.data());
    int const file_status = MPI_File_open(_mpi_comm, filename_char,
            MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if(file_status != 0)
    {
        ERR("Error opening file %s. MPI error code %d", filename.c_str(), file_status);
        return false;
    }

    // Read data
    char file_mode[] = "native";
    MPI_File_set_view(file, offset, type, type, file_mode, MPI_INFO_NULL);
    MPI_File_read(file, data.data(), data.size(), type, MPI_STATUS_IGNORE);
    MPI_File_close(&file);

    return true;
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::readBinary(
        const std::string &file_name_base)
{
    //----------------------------------------------------------------------------------
    // Read headers

    const std::string fname_header = file_name_base +  "_partitioned_msh_";
    const std::string fname_num_p_ext = std::to_string(_mpi_comm_size) + ".bin";

    // Array that contains integer numbers of the partition data header.
    // See description of readBinary() function.
    std::array<long, 14> mesh_controls;

    if (!readBinaryDataFromFile(fname_header + "cfg" + fname_num_p_ext,
            static_cast<MPI_Offset>(
                static_cast<unsigned>(_mpi_rank) * mesh_controls.size() * sizeof(long)),
            MPI_LONG, mesh_controls))
        return nullptr;

    _num_nodes_part = mesh_controls[0];
    _num_regular_elems_part = mesh_controls[2];
    _num_ghost_elems_part = mesh_controls[3];

    //----------------------------------------------------------------------------------
    // Read Nodes
    std::vector<NodeData> nodes(static_cast<std::size_t>(_num_nodes_part));

    if (!readBinaryDataFromFile(fname_header + "nod" + fname_num_p_ext,
             static_cast<MPI_Offset>(mesh_controls[10]), _mpi_node_type, nodes))
        return nullptr;

    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<std::size_t> glb_node_ids;
    setNodes(nodes, mesh_nodes, glb_node_ids);

    //----------------------------------------------------------------------------------
    // Read non-ghost elements

    std::vector<long> elem_data(static_cast<std::size_t>(
        _num_regular_elems_part + mesh_controls[8]));
    if (!readBinaryDataFromFile(fname_header +"ele" + fname_num_p_ext,
            static_cast<MPI_Offset>(mesh_controls[11]), MPI_LONG, elem_data))
        return nullptr;

    std::vector<MeshLib::Element*> mesh_elems(_num_regular_elems_part + _num_ghost_elems_part);
    setElements(mesh_nodes, elem_data, mesh_elems);

    //----------------------------------------------------------------------------------
    //Read ghost element
    std::vector<long> ghost_elem_data(static_cast<std::size_t>(
        _num_ghost_elems_part + mesh_controls[9]));

    if (!readBinaryDataFromFile(fname_header + "ele_g" + fname_num_p_ext,
            static_cast<MPI_Offset>(mesh_controls[12]), MPI_LONG, ghost_elem_data))
        return nullptr;

    const bool process_ghost = true;
    setElements(mesh_nodes, ghost_elem_data, mesh_elems, process_ghost);

    //----------------------------------------------------------------------------------
    //Create a mesh and return it
    const std::size_t n_nghost_elem = static_cast<std::size_t>(mesh_controls[2]);
    const unsigned n_global_base_nodes = static_cast<unsigned>(mesh_controls[6]);
    const unsigned n_global_nodes = static_cast<unsigned>(mesh_controls[7]);
    const unsigned n_base_nodes = static_cast<unsigned>(mesh_controls[1]);
    const unsigned n_active_base_nodes = static_cast<unsigned>(mesh_controls[4]);
    const unsigned n_active_nodes = static_cast<unsigned>(mesh_controls[5]);

    return new MeshLib::NodePartitionedMesh(
            BaseLib::extractBaseName(file_name_base) + std::to_string(_mpi_comm_size),
            mesh_nodes, glb_node_ids,
            mesh_elems, n_nghost_elem,
            n_global_base_nodes, n_global_nodes,
            n_base_nodes, n_active_base_nodes, n_active_nodes);
}

bool NodePartitionedMeshReader::openASCIIFiles(std::string const& file_name_base,
        std::ifstream& is_cfg, std::ifstream& is_node, std::ifstream& is_elem) const
{
    const std::string fname_header = file_name_base +  "_partitioned_";
    const std::string fname_num_p_ext = std::to_string(_mpi_comm_size) + ".msh";

    {   // Configuration.
        std::string const filename = fname_header + "cfg"+ fname_num_p_ext;
        is_cfg.open(filename);

        if( !is_cfg.good() )
        {
            ERR("Error opening file %s for input.", filename.c_str());
            return false;
        }

        std::string tmp_line;
        std::getline(is_cfg, tmp_line);
        int num_parts;
        is_cfg >> num_parts >> std::ws;

        if(num_parts != _mpi_comm_size)
        {
            ERR("Aborting computation because of number of cores"
                "/ subdomains mismatch.");
            return false;
        }
    }

    {   // Nodes.
        std::string const filename = fname_header + "nodes_" + fname_num_p_ext;
        is_node.open(filename);
        if( !is_node.good() )
        {
            ERR("Error opening file %s for input.", filename.c_str());
            return false;
        }
    }

    {   // Elements.
        std::string const filename = fname_header + "elems_" + fname_num_p_ext;
        is_elem.open(filename);
        if( !is_elem.good() )
        {
            ERR("Error opening file %s for input.", filename.c_str());
            return false;
        }
    }

    return true;
}

void NodePartitionedMeshReader::readCastNodesASCII(std::ifstream& is_node,
        const int part_id, std::vector<MeshLib::Node*> &mesh_nodes,
        std::vector<std::size_t> &glb_node_ids) const
{
    int const message_tag = 0;

    //----------------------------------------------------------------------------------
    // Read Nodes
    std::vector<NodeData> nodes(static_cast<std::size_t>(_num_nodes_part));

    if(_mpi_rank == 0)
    {
        for(std::size_t k=0; k<_num_nodes_part; k++)
        {
            NodeData &node = nodes[k];
            is_node >> node.index
                    >> node.x >> node.y >> node.z >> std::ws;
        }

        if(part_id == 0)
        {
            setNodes(nodes, mesh_nodes, glb_node_ids);
        }
        else
        {
            MPI_Send(nodes.data(), _num_nodes_part, _mpi_node_type, part_id, message_tag, _mpi_comm);
        }
    }
    else if(part_id > 0 && _mpi_rank == part_id)
    {
        MPI_Recv(nodes.data(), _num_nodes_part, _mpi_node_type, 0, message_tag, _mpi_comm, MPI_STATUS_IGNORE);
        setNodes(nodes, mesh_nodes, glb_node_ids);
    }
}

void NodePartitionedMeshReader::readCastElemsASCII(std::ifstream& is_elem,
        const int part_id, const int data_size, const bool process_ghost,
        const std::vector<MeshLib::Node*> &mesh_nodes,
        std::vector<MeshLib::Element*> &mesh_elems) const
{
    int const message_tag = 0;

    std::vector<long> elem_data(static_cast<std::size_t>(data_size));
    if(_mpi_rank == 0)
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
            MPI_Send(elem_data.data(), data_size, MPI_LONG, part_id, message_tag, _mpi_comm);
        }
    }
    else if(part_id > 0 && _mpi_rank == part_id)
    {
        MPI_Recv(elem_data.data(), data_size, MPI_LONG, 0, message_tag, _mpi_comm, MPI_STATUS_IGNORE);

        if(!process_ghost)
            mesh_elems.resize(_num_regular_elems_part + _num_ghost_elems_part);
        setElements(mesh_nodes, elem_data, mesh_elems, process_ghost);
    }
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::readASCII(
        const std::string &file_name_base)
{
    std::ifstream is_cfg;
    std::ifstream is_node;
    std::ifstream is_elem;

    bool file_opened = false;
    if(_mpi_rank == 0)
    {
        file_opened = openASCIIFiles(file_name_base, is_cfg, is_node, is_elem);

        if(!file_opened)
        {
            is_cfg.close();
            is_node.close();
            is_elem.close();
        }
    }

    MPI_Bcast(&file_opened, 1, MPI_INT, 0, _mpi_comm);

    if(!file_opened)
        return nullptr;

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
    std::array<long, 11> mesh_controls;

    MeshLib::NodePartitionedMesh *np_mesh = nullptr;
    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<std::size_t> glb_node_ids;
    std::vector<MeshLib::Element*> mesh_elems;

    for(int i = 0; i < _mpi_comm_size; i++)
    {
        if(_mpi_rank == 0)
        {
            INFO("-->Parallel reading the partitioned mesh: ");

            for(long& v : mesh_controls)
                is_cfg >> v;
            is_cfg >> std::ws;
        }

        MPI_Bcast(mesh_controls.data(), mesh_controls.size(), MPI_LONG, 0, _mpi_comm);

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

        if(_mpi_rank == i)
        {
            //----------------------------------------------------------------------------------
            //Create a mesh
            const std::size_t n_nghost_elem = static_cast<std::size_t>(mesh_controls[2]);
            const unsigned n_global_base_nodes = static_cast<unsigned>(mesh_controls[6]);
            const unsigned n_global_nodes = static_cast<unsigned>(mesh_controls[7]);
            const unsigned n_base_nodes = static_cast<unsigned>(mesh_controls[1]);
            const unsigned n_active_base_nodes = static_cast<unsigned>(mesh_controls[4]);
            const unsigned n_active_nodes = static_cast<unsigned>(mesh_controls[5]);

            np_mesh = new  MeshLib::NodePartitionedMesh(
                    BaseLib::extractBaseName(file_name_base) + std::to_string(_mpi_comm_size),
                    mesh_nodes, glb_node_ids,
                    mesh_elems, n_nghost_elem,
                    n_global_base_nodes, n_global_nodes,
                    n_base_nodes, n_active_base_nodes, n_active_nodes);
        }
    }

    if(_mpi_rank == 0)
    {
        is_cfg.close();
        is_node.close();
        is_elem.close();
    }

    MPI_Barrier(_mpi_comm);

    return  np_mesh;
}

void NodePartitionedMeshReader::readElementASCII(std::ifstream &ins,
        std::vector<long>& elem_data, const bool ghost) const
{
    // Set number of elements.
    const long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    long counter = ne;
    for(long j=0; j<ne; j++)
    {
        elem_data[j] = counter;
        ins >> elem_data[counter++];  //mat. idx
        ins >> elem_data[counter++];  //type
        ins >> elem_data[counter];  //nnodes
        const long nn_e =  elem_data[counter++];
        for(long k=0; k<nn_e; k++)
            ins >> elem_data[counter++];
    }
}

void NodePartitionedMeshReader::setNodes(const std::vector<NodeData> &node_data,
        std::vector<MeshLib::Node*> &mesh_node,
        std::vector<std::size_t> &glb_node_ids) const
{
    mesh_node.resize( _num_nodes_part );
    glb_node_ids.resize( _num_nodes_part );

    for(std::size_t i=0; i<mesh_node.size(); i++)
    {
        NodeData const& nd = node_data[i];
        glb_node_ids[i] = nd.index;
        mesh_node[i] = new MeshLib::Node(nd.x, nd.y, nd.z, i);
    }
}

void NodePartitionedMeshReader::setElements(
        const std::vector<MeshLib::Node*> &mesh_nodes,
        const std::vector<long> &elem_data,
        std::vector<MeshLib::Element*> &mesh_elems, const bool ghost) const
{
    // Number of elements, ether ghost or regular
    const long ne = ghost ? _num_ghost_elems_part : _num_regular_elems_part;
    const long id_offset = ghost ? _num_regular_elems_part : 0;

    for(std::size_t i=0; i < ne; i++)
    {
        std::size_t counter = elem_data[i];

        const unsigned mat_idx = static_cast<unsigned>( elem_data[counter++] );
        const long e_type = elem_data[counter++];
        const long nnodes = elem_data[counter++];

        MeshLib::Node **elem_nodes = new MeshLib::Node*[nnodes];
        for(std::size_t k=0; k<nnodes; k++)
            elem_nodes[k] = mesh_nodes[ elem_data[counter++] ];

        MeshLib::Element *elem = nullptr;

        // The element types below are defined by the mesh_partition tool
        // available at https://github.com/ufz/mesh_partition .
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

void NodePartitionedMeshReader::registerNodeDataMpiType()
{
    int const count = 2;
    int blocks[count] = {1, 3};
    MPI_Datatype types[count] = {MPI_UNSIGNED_LONG, MPI_DOUBLE};
    MPI_Aint displacements[count] = {0, sizeof(NodeData::index)};

    MPI_Type_create_struct(count, blocks, displacements, types, &_mpi_node_type);
    MPI_Type_commit(&_mpi_node_type);
}

}   // namespace FileIO
