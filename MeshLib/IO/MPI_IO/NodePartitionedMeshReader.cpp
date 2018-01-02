/*!
  \file NodePartitionedMeshReader.cpp
  \author Wenqing Wang
  \date   2014.08
  \brief  Define members of class NodePartitionedMeshReader to read node-wise
  partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

 */

#include "NodePartitionedMeshReader.h"

#include <logog/include/logog.hpp>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Properties.h"

// Check if the value can by converted to given type without overflow.
template <typename VALUE, typename TYPE>
bool
is_safely_convertable(VALUE const& value)
{
    bool const result = value <= std::numeric_limits<TYPE>::max();
    if (!result)
    {
        ERR("The value %d is too large for conversion.", value);
        ERR("Maximum available size is %d.", std::numeric_limits<TYPE>::max());
    }
    return result;
}

namespace MeshLib
{
namespace IO
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

    MeshLib::NodePartitionedMesh* mesh = nullptr;

    // Always try binary file first
    std::string const fname_new = file_name_base + "_partitioned_msh_cfg" +
        std::to_string(_mpi_comm_size) + ".bin";

    if(!BaseLib::IsFileExisting(fname_new)) // doesn't exist binary file.
    {
        INFO("Reading ASCII mesh file ...");

        mesh = readASCII(file_name_base);
    }
    else
    {
        INFO("Reading binary mesh file ...");

        mesh = readBinary(file_name_base);
    }

    INFO("[time] Reading the mesh took %f s.", timer.elapsed());

    MPI_Barrier(_mpi_comm);

    return mesh;
}

template <typename DATA>
bool
NodePartitionedMeshReader::readBinaryDataFromFile(std::string const& filename,
    MPI_Offset offset, MPI_Datatype type, DATA& data) const
{
    // Check container size
    if (!is_safely_convertable<std::size_t, int>(data.size()))
    {
        ERR("The container size is too large for MPI_File_read() call.");
        return false;
    }

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
    // The static cast is checked above.
    MPI_File_read(file, data.data(), static_cast<int>(data.size()), type,
        MPI_STATUS_IGNORE);
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
    std::vector<NodeData> nodes(_mesh_info.nodes);

    if (!readBinaryDataFromFile(fname_header + "nod" + fname_num_p_ext,
        static_cast<MPI_Offset>(_mesh_info.offset[2]), _mpi_node_type, nodes))
        return nullptr;

    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<unsigned long> glb_node_ids;
    setNodes(nodes, mesh_nodes, glb_node_ids);

    //----------------------------------------------------------------------------------
    // Read non-ghost elements
    std::vector<unsigned long> elem_data(
        _mesh_info.regular_elements + _mesh_info.offset[0]);
    if (!readBinaryDataFromFile(fname_header + "ele" + fname_num_p_ext,
        static_cast<MPI_Offset>(_mesh_info.offset[3]), MPI_LONG, elem_data))
        return nullptr;

    std::vector<MeshLib::Element*> mesh_elems(
        _mesh_info.regular_elements + _mesh_info.ghost_elements);
    setElements(mesh_nodes, elem_data, mesh_elems);

    //----------------------------------------------------------------------------------
    //Read ghost element
    std::vector<unsigned long> ghost_elem_data(
        _mesh_info.ghost_elements + _mesh_info.offset[1]);

    if (!readBinaryDataFromFile(fname_header + "ele_g" + fname_num_p_ext,
        static_cast<MPI_Offset>(_mesh_info.offset[4]), MPI_LONG, ghost_elem_data))
        return nullptr;

    const bool process_ghost = true;
    setElements(mesh_nodes, ghost_elem_data, mesh_elems, process_ghost);

    //----------------------------------------------------------------------------------
    // read the properties
    MeshLib::Properties p(readPropertiesBinary(file_name_base));

    return newMesh(BaseLib::extractBaseName(file_name_base), mesh_nodes,
                   glb_node_ids, mesh_elems, p);
}

MeshLib::Properties NodePartitionedMeshReader::readPropertiesBinary(
    const std::string& file_name_base) const
{
    MeshLib::Properties p;
    readPropertiesBinary(file_name_base, MeshLib::MeshItemType::Node, p);
    readPropertiesBinary(file_name_base, MeshLib::MeshItemType::Cell, p);
    return p;
}

void NodePartitionedMeshReader::readPropertiesBinary(
    const std::string& file_name_base, MeshLib::MeshItemType t,
    MeshLib::Properties& p) const
{
    std::string const item_type =
        t == MeshLib::MeshItemType::Node ? "node" : "cell";
    const std::string fname_cfg = file_name_base + "_partitioned_" + item_type +
                                  "_properties_cfg" +
                                  std::to_string(_mpi_comm_size) + ".bin";
    std::ifstream is(fname_cfg.c_str(), std::ios::binary | std::ios::in);
    if (!is)
    {
        WARN("Could not open file '%s'.\n"
             "\tYou can ignore this warning if the mesh does not contain %s-"
             "wise property data.", fname_cfg.c_str(), item_type.data());
        return;
    }
    std::size_t number_of_properties = 0;
    is.read(reinterpret_cast<char*>(&number_of_properties), sizeof(std::size_t));
    std::vector<boost::optional<MeshLib::IO::PropertyVectorMetaData>> vec_pvmd(
        number_of_properties);
    for (std::size_t i(0); i < number_of_properties; ++i)
    {
        vec_pvmd[i] = MeshLib::IO::readPropertyVectorMetaData(is);
        if (!vec_pvmd[i])
        {
            OGS_FATAL(
                "Error in NodePartitionedMeshReader::readPropertiesBinary: "
                "Could not read the meta data for the PropertyVector %d",
                i);
        }
    }
    for (std::size_t i(0); i < number_of_properties; ++i)
    {
        DBUG("[%d] +++++++++++++", _mpi_rank);
        MeshLib::IO::writePropertyVectorMetaData(*(vec_pvmd[i]));
        DBUG("[%d] +++++++++++++", _mpi_rank);
    }
    auto pos = is.tellg();
    auto offset =
        static_cast<long>(pos) +
        static_cast<long>(_mpi_rank *
                          sizeof(MeshLib::IO::PropertyVectorPartitionMetaData));
    is.seekg(offset);
    boost::optional<MeshLib::IO::PropertyVectorPartitionMetaData> pvpmd(
        MeshLib::IO::readPropertyVectorPartitionMetaData(is));
    bool pvpmd_read_ok = static_cast<bool>(pvpmd);
    bool all_pvpmd_read_ok;
    MPI_Allreduce(&pvpmd_read_ok, &all_pvpmd_read_ok, 1, MPI_C_BOOL, MPI_LOR,
                  _mpi_comm);
    if (!all_pvpmd_read_ok)
    {
        OGS_FATAL(
            "Error in NodePartitionedMeshReader::readPropertiesBinary: "
            "Could not read the partition meta data for the mpi process %d",
            _mpi_rank);
    }
    DBUG("[%d] offset in the PropertyVector: %d", _mpi_rank, pvpmd->offset);
    DBUG("[%d] %d tuples in partition.", _mpi_rank, pvpmd->number_of_tuples);
    is.close();

    const std::string fname_val = file_name_base + "_partitioned_" + item_type +
                                  "_properties_val" +
                                  std::to_string(_mpi_comm_size) + ".bin";
    is.open(fname_val.c_str(), std::ios::binary | std::ios::in);
    if (!is)
    {
        ERR("Could not open file '%s'\n."
            "\tYou can ignore this warning if the mesh does not contain %s-"
             "wise property data.", fname_val.c_str(), item_type.data());
    }

    readDomainSpecificPartOfPropertyVectors(vec_pvmd, *pvpmd, t, is, p);
}

void NodePartitionedMeshReader::readDomainSpecificPartOfPropertyVectors(
    std::vector<boost::optional<MeshLib::IO::PropertyVectorMetaData>> const&
        vec_pvmd,
    MeshLib::IO::PropertyVectorPartitionMetaData const& pvpmd,
    MeshLib::MeshItemType t,
    std::istream& is,
    MeshLib::Properties& p) const
{
    unsigned long global_offset = 0;
    std::size_t const number_of_properties = vec_pvmd.size();
    for (std::size_t i(0); i < number_of_properties; ++i)
    {
        DBUG("[%d] global offset: %d, offset within the PropertyVector: %d.",
             _mpi_rank, global_offset,
             global_offset +
                 pvpmd.offset * vec_pvmd[i]->data_type_size_in_bytes);
        if (vec_pvmd[i]->is_int_type)
        {
            if (vec_pvmd[i]->is_data_type_signed)
            {
                if (vec_pvmd[i]->data_type_size_in_bytes == sizeof(int))
                    createPropertyVectorPart<int>(is, *vec_pvmd[i], pvpmd, t,
                                                  global_offset, p);
                if (vec_pvmd[i]->data_type_size_in_bytes == sizeof(long))
                    createPropertyVectorPart<long>(is, *vec_pvmd[i], pvpmd, t,
                                                   global_offset, p);
            }
            else
            {
                if (vec_pvmd[i]->data_type_size_in_bytes ==
                    sizeof(unsigned int))
                    createPropertyVectorPart<unsigned int>(
                        is, *vec_pvmd[i], pvpmd, t, global_offset, p);
                if (vec_pvmd[i]->data_type_size_in_bytes ==
                    sizeof(unsigned long))
                    createPropertyVectorPart<unsigned long>(
                        is, *vec_pvmd[i], pvpmd, t, global_offset, p);
            }
        }
        else
        {
            if (vec_pvmd[i]->data_type_size_in_bytes == sizeof(float))
                createPropertyVectorPart<float>(is, *vec_pvmd[i], pvpmd, t,
                                                global_offset, p);
            if (vec_pvmd[i]->data_type_size_in_bytes == sizeof(double))
                createPropertyVectorPart<double>(is, *vec_pvmd[i], pvpmd, t,
                                                 global_offset, p);
        }
        global_offset += vec_pvmd[i]->data_type_size_in_bytes *
                         vec_pvmd[i]->number_of_tuples *
                         vec_pvmd[i]->number_of_components;
    }
}

bool NodePartitionedMeshReader::openASCIIFiles(
    std::string const& file_name_base, std::ifstream& is_cfg,
    std::ifstream& is_node, std::ifstream& is_elem) const
{
    const std::string fname_header = file_name_base +  "_partitioned_";
    const std::string fname_num_p_ext = std::to_string(_mpi_comm_size) + ".msh";

    {   // Configuration.
        std::string const filename = fname_header + "cfg" + fname_num_p_ext;
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

bool NodePartitionedMeshReader::readCastNodesASCII(std::ifstream& is_node,
    const int part_id, std::vector<MeshLib::Node*> &mesh_nodes,
    std::vector<unsigned long> &glb_node_ids) const
{
    int const message_tag = 0;

    // MPI_Send/Recv can handle only int. Check for overflow.
    if (!is_safely_convertable<unsigned long, int>(_mesh_info.nodes))
    {
        ERR("Too large number of nodes to read.");
        return false;
    }

    std::vector<NodeData> nodes(_mesh_info.nodes);

    if(_mpi_rank == 0)
    {
        for(unsigned long k = 0; k < _mesh_info.nodes; k++)
        {
            NodeData &node = nodes[k];
            is_node >> node.index >> node.x >> node.y >> node.z >> std::ws;
        }

        if(part_id == 0)
            setNodes(nodes, mesh_nodes, glb_node_ids);
        else
            MPI_Send(nodes.data(), static_cast<int>(_mesh_info.nodes),
                _mpi_node_type, part_id, message_tag, _mpi_comm);
    }
    else if(_mpi_rank == part_id)
    {
        MPI_Recv(nodes.data(), static_cast<int>(_mesh_info.nodes),
            _mpi_node_type, 0, message_tag, _mpi_comm, MPI_STATUS_IGNORE);
        setNodes(nodes, mesh_nodes, glb_node_ids);
    }

    return true;
}

bool NodePartitionedMeshReader::readCastElemsASCII(std::ifstream& is_elem,
    const int part_id, const std::size_t data_size, const bool process_ghost,
    const std::vector<MeshLib::Node*> &mesh_nodes,
    std::vector<MeshLib::Element*> &mesh_elems) const
{
    int const message_tag = 0;

    // MPI_Send/Recv can handle only int. Check for overflow.
    if (!is_safely_convertable<std::size_t, int>(data_size))
    {
        ERR("Too large number of elements to read.");
        return false;
    }

    std::vector<unsigned long> elem_data(data_size);
    if(_mpi_rank == 0)
    {
        readElementASCII(is_elem, elem_data, process_ghost);

        if(part_id == 0)
        {
            if(!process_ghost)
                mesh_elems.resize(
                    _mesh_info.regular_elements + _mesh_info.ghost_elements);
            setElements(mesh_nodes, elem_data, mesh_elems, process_ghost);
        }
        else
            MPI_Send(elem_data.data(), static_cast<int>(data_size), MPI_LONG,
                part_id, message_tag, _mpi_comm);
    }
    else if(_mpi_rank == part_id)
    {
        MPI_Recv(elem_data.data(), static_cast<int>(data_size), MPI_LONG,
            0, message_tag, _mpi_comm, MPI_STATUS_IGNORE);

        if(!process_ghost)
            mesh_elems.resize(
                _mesh_info.regular_elements + _mesh_info.ghost_elements);
        setElements(mesh_nodes, elem_data, mesh_elems, process_ghost);
    }

    return true;
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

    MeshLib::NodePartitionedMesh* np_mesh = nullptr;
    std::vector<MeshLib::Node*> mesh_nodes;
    std::vector<unsigned long> glb_node_ids;
    std::vector<MeshLib::Element*> mesh_elems;

    for(int i = 0; i < _mpi_comm_size; i++)
    {
        if(_mpi_rank == 0)
        {
            // Read first part into _mesh_info which is equal with the binary
            // structure.
            for(unsigned long j = 0; j < 10; ++j)
                is_cfg >> *(_mesh_info.data() + j);
            // The last positon is the extra_flag.
            is_cfg >> _mesh_info.extra_flag;
            is_cfg >> std::ws;
        }

        MPI_Bcast(_mesh_info.data(), static_cast<int>(_mesh_info.size()),
            MPI_LONG, 0, _mpi_comm);

        //---------------------------------------------------------------------
        // Read Nodes
        if (!readCastNodesASCII(is_node, i, mesh_nodes, glb_node_ids))
            break;

        //---------------------------------------------------------------------
        // Read elements
        if (!readCastElemsASCII(is_elem, i,
            _mesh_info.regular_elements + _mesh_info.offset[0],
            false, mesh_nodes, mesh_elems))
            break;

        //---------------------------------------------------------------------
        // Ghost elements
        if (!readCastElemsASCII(is_elem, i,
            _mesh_info.ghost_elements + _mesh_info.offset[1],
            true, mesh_nodes, mesh_elems))
            break;

        if(_mpi_rank == i) {
            // reading ascii properties is not implemented
            MeshLib::Properties properties;
            np_mesh = newMesh(BaseLib::extractBaseName(file_name_base),
                    mesh_nodes, glb_node_ids, mesh_elems, properties);
        }
    }

    if(_mpi_rank == 0)
    {
        is_cfg.close();
        is_node.close();
        is_elem.close();
    }

    MPI_Barrier(_mpi_comm);

    return np_mesh;
}

MeshLib::NodePartitionedMesh* NodePartitionedMeshReader::newMesh(
    std::string const& mesh_name,
    std::vector<MeshLib::Node*> const& mesh_nodes,
    std::vector<unsigned long> const& glb_node_ids,
    std::vector<MeshLib::Element*> const& mesh_elems,
    MeshLib::Properties const& properties) const
{
    return new MeshLib::NodePartitionedMesh(
        mesh_name + std::to_string(_mpi_comm_size),
        mesh_nodes, glb_node_ids, mesh_elems,
        properties,
        _mesh_info.global_base_nodes,
        _mesh_info.global_nodes,
        _mesh_info.base_nodes,
        _mesh_info.active_base_nodes,
        _mesh_info.active_nodes);
}

void NodePartitionedMeshReader::readElementASCII(std::ifstream &ins,
    std::vector<unsigned long>& elem_data, const bool ghost) const
{
    // Set number of elements.
    unsigned long const ne =
        ghost ? _mesh_info.ghost_elements : _mesh_info.regular_elements;
    unsigned long id_offset_elem = ne;
    for(unsigned long j = 0; j < ne; j++)
    {
        elem_data[j] = id_offset_elem;
        ins >> elem_data[id_offset_elem++];  //mat. idx
        ins >> elem_data[id_offset_elem++];  //type
        ins >> elem_data[id_offset_elem];  //nnodes
        unsigned long const nn_e =  elem_data[id_offset_elem++];
        for(unsigned long k = 0; k < nn_e; k++)
            ins >> elem_data[id_offset_elem++];
    }
}

void NodePartitionedMeshReader::setNodes(const std::vector<NodeData> &node_data,
    std::vector<MeshLib::Node*> &mesh_node,
    std::vector<unsigned long> &glb_node_ids) const
{
    mesh_node.resize(_mesh_info.nodes);
    glb_node_ids.resize(_mesh_info.nodes);

    for(std::size_t i = 0; i < mesh_node.size(); i++)
    {
        NodeData const& nd = node_data[i];
        glb_node_ids[i] = nd.index;
        mesh_node[i] = new MeshLib::Node(nd.x, nd.y, nd.z, i);
    }
}

void NodePartitionedMeshReader::setElements(
    const std::vector<MeshLib::Node*> &mesh_nodes,
    const std::vector<unsigned long> &elem_data,
    std::vector<MeshLib::Element*> &mesh_elems, const bool ghost) const
{
    // Number of elements, ether ghost or regular
    unsigned long const ne =
        ghost ? _mesh_info.ghost_elements : _mesh_info.regular_elements;
    unsigned long const id_offset_ghost =
        ghost ? _mesh_info.regular_elements : 0;

    for(unsigned long i = 0; i < ne; i++)
    {
        unsigned long id_offset_elem = elem_data[i];

        const unsigned mat_idx = static_cast<unsigned>( elem_data[id_offset_elem++] );
        const unsigned long e_type = elem_data[id_offset_elem++];
        unsigned long const nnodes = elem_data[id_offset_elem++];

        MeshLib::Node** elem_nodes = new MeshLib::Node*[nnodes];
        for(unsigned long k = 0; k < nnodes; k++)
            elem_nodes[k] = mesh_nodes[ elem_data[id_offset_elem++] ];

        // The element types below are defined by the mesh_partition tool
        // available at https://github.com/ufz/mesh_partition .
        switch(e_type)
        {
        case 2:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Line(elem_nodes, mat_idx);
            break;
        case 6:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Quad(elem_nodes, mat_idx);
            break;
        case 11:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Hex(elem_nodes, mat_idx);
            break;
        case 4:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Tri(elem_nodes, mat_idx);
            break;
        case 9:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Tet(elem_nodes, mat_idx);
            break;
        case 14:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Prism(elem_nodes, mat_idx);
            break;
        case 17:
            mesh_elems[i + id_offset_ghost] = new MeshLib::Pyramid(elem_nodes, mat_idx);
            break;
        }
    }
}
}   // namespace IO
}   // namespace MeshLib
