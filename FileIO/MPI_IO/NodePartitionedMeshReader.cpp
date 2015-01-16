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

// Check if the value can by converted to given type without overflow.
template <typename VALUE, typename TYPE>
bool
is_safely_convertable(VALUE const& value)
{
    bool const result = value <= std::numeric_limits<TYPE>::max();
    if (!result) {
        ERR("The value %d is too large for conversion.", value);
        ERR("Maximum available size is %d.", std::numeric_limits<TYPE>::max());
    }
    return result;
}

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


void NodePartitionedMeshReader::registerNodeDataMpiType()
{
    int const count = 2;
    int blocks[count] = {1, 3};
    MPI_Datatype types[count] = {MPI_UNSIGNED_LONG, MPI_DOUBLE};
    MPI_Aint displacements[count] = {0, sizeof(NodeData::index)};

    MPI_Type_create_struct(count, blocks, displacements, types, &_mpi_node_type);
    MPI_Type_commit(&_mpi_node_type);
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

    if (!readBinaryDataFromFile(fname_header + "cfg" + fname_num_p_ext,
            static_cast<MPI_Offset>(
                static_cast<unsigned>(_mpi_rank) * sizeof(_mesh_info)),
            MPI_LONG, _mesh_info))
        return nullptr;

    //----------------------------------------------------------------------------------
    // Read Nodes
    std::vector<NodeData> nodes(static_cast<std::size_t>(_mesh_info.nodes));

    if (!readBinaryDataFromFile(fname_header + "nod" + fname_num_p_ext,
             static_cast<MPI_Offset>(_mesh_info.offset[2]), _mpi_node_type, nodes))
        return nullptr;

    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<std::size_t> glb_node_ids;
    setNodes(nodes, mesh_nodes, glb_node_ids);

    //----------------------------------------------------------------------------------
    // Read non-ghost elements

    std::vector<long> elem_data(static_cast<std::size_t>(
        _mesh_info.regular_elements + _mesh_info.offset[0]));
    if (!readBinaryDataFromFile(fname_header +"ele" + fname_num_p_ext,
            static_cast<MPI_Offset>(_mesh_info.offset[3]), MPI_LONG, elem_data))
        return nullptr;

    std::vector<MeshLib::Element*> mesh_elems(
            _mesh_info.regular_elements + _mesh_info.ghost_elements);
    setElements(mesh_nodes, elem_data, mesh_elems);

    //----------------------------------------------------------------------------------
    //Read ghost element
    std::vector<long> ghost_elem_data(static_cast<std::size_t>(
        _mesh_info.ghost_elements + _mesh_info.offset[1]));

    if (!readBinaryDataFromFile(fname_header + "ele_g" + fname_num_p_ext,
            static_cast<MPI_Offset>(_mesh_info.offset[4]), MPI_LONG, ghost_elem_data))
        return nullptr;

    const bool process_ghost = true;
    setElements(mesh_nodes, ghost_elem_data, mesh_elems, process_ghost);

    //----------------------------------------------------------------------------------
    return newMesh(BaseLib::extractBaseName(file_name_base),
            mesh_nodes, glb_node_ids, mesh_elems);
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
    std::vector<NodeData> nodes(static_cast<std::size_t>(_mesh_info.nodes));

    if(_mpi_rank == 0)
    {
        for(std::size_t k=0; k < static_cast<std::size_t>(_mesh_info.nodes); k++)
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
            MPI_Send(nodes.data(), _mesh_info.nodes, _mpi_node_type, part_id,
                    message_tag, _mpi_comm);
        }
    }
    else if(part_id > 0 && _mpi_rank == part_id)
    {
        MPI_Recv(nodes.data(), _mesh_info.nodes, _mpi_node_type, 0,
                message_tag, _mpi_comm, MPI_STATUS_IGNORE);
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
                mesh_elems.resize(_mesh_info.regular_elements +
                                    _mesh_info.ghost_elements);
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
            mesh_elems.resize(_mesh_info.regular_elements +
                                _mesh_info.ghost_elements);
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

    MeshLib::NodePartitionedMesh *np_mesh = nullptr;
    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<std::size_t> glb_node_ids;
    std::vector<MeshLib::Element*> mesh_elems;

    for(int i = 0; i < _mpi_comm_size; i++)
    {
        if(_mpi_rank == 0)
        {
            INFO("-->Parallel reading the partitioned mesh: ");

            // Read first part into _mesh_info which is equal with the binary
            // structure.
            for(std::size_t j = 0; j < 10; ++j)
                is_cfg >> *(_mesh_info.data() + j);
            // The last positon is the extra_flag.
            is_cfg >> _mesh_info.extra_flag;
            is_cfg >> std::ws;
        }

        MPI_Bcast(_mesh_info.data(), _mesh_info.size(), MPI_LONG, 0, _mpi_comm);

        //----------------------------------------------------------------------------------
        // Read Nodes
        readCastNodesASCII(is_node, i, mesh_nodes, glb_node_ids);

        //----------------------------------------------------------------------------------
        // Read elements
        readCastElemsASCII(is_elem, i,
            _mesh_info.regular_elements + _mesh_info.offset[0],
            false, mesh_nodes, mesh_elems);

        //-------------------------------------------------------------------------
        // Ghost elements
        readCastElemsASCII(is_elem, i,
            _mesh_info.ghost_elements + _mesh_info.offset[1],
            true, mesh_nodes, mesh_elems);

        if(_mpi_rank == i)
        {
            np_mesh = newMesh(BaseLib::extractBaseName(file_name_base),
                    mesh_nodes, glb_node_ids, mesh_elems);
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

MeshLib::NodePartitionedMesh*
NodePartitionedMeshReader::newMesh(
    std::string const& mesh_name,
    std::vector<MeshLib::Node*> const& mesh_nodes,
    std::vector<std::size_t> const& glb_node_ids,
    std::vector<MeshLib::Element*> const& mesh_elems) const
{
    return new MeshLib::NodePartitionedMesh(
        mesh_name + std::to_string(_mpi_comm_size),
        mesh_nodes, glb_node_ids, mesh_elems,
        static_cast<std::size_t>(_mesh_info.regular_elements),
        static_cast<unsigned>(_mesh_info.global_base_nodes),
        static_cast<unsigned>(_mesh_info.global_nodes),
        static_cast<unsigned>(_mesh_info.base_nodes),
        static_cast<unsigned>(_mesh_info.active_base_nodes),
        static_cast<unsigned>(_mesh_info.active_nodes));
}

void NodePartitionedMeshReader::readElementASCII(std::ifstream &ins,
        std::vector<long>& elem_data, const bool ghost) const
{
    // Set number of elements.
    const long ne = ghost ? _mesh_info.ghost_elements : _mesh_info.regular_elements;
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
    mesh_node.resize(_mesh_info.nodes);
    glb_node_ids.resize(_mesh_info.nodes);

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
    const long ne = ghost ? _mesh_info.ghost_elements : _mesh_info.regular_elements;
    const long id_offset = ghost ? _mesh_info.regular_elements : 0;

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

}   // namespace FileIO
