/*!
  \file MeshPartitioning.cpp
  \date   2016.05

  \brief  Define the members of class MeshPartitioning

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "MeshPartitioning.h"

#include <iomanip>
#include <stdio.h> // for binary output

#include "logog/include/logog.hpp"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
void MeshPartitioning :: write2METIS(const std::string& file_name)
{
    std::ofstream os(file_name.data(), std::ios::trunc);
    os << _elements.size() <<" \n";
    for (const auto elem : _elements)
    {
        for (unsigned j=0; j<elem->getNNodes(); j++)
        {
            os << elem->getNodeIndex(j) + 1 <<" ";
        }
        os << std::endl;
    }
}

void MeshPartitioning :: partitionByNodeMETIS(const std::string& file_name_base,
        const unsigned npartitions,
        const bool output_binary)
{
    const std::string npartitions_str = std::to_string(npartitions);

    // -- Read partitioned mesh data from METIS
    std::string fname_parts = file_name_base + ".mesh.npart." + npartitions_str;

    std::ifstream npart_in(fname_parts.data());
    if (!npart_in.is_open())
    {
        ERR("Error: cannot open file %s. It may not exist! ", fname_parts.data());
        abort();
    }

    const std::size_t nnodes = _nodes.size();

    std::vector<std::size_t> nodes_partition_ids;
    nodes_partition_ids.reserve(nnodes);

    for (std::size_t i=0; i<nnodes; i++)
    {
        std::size_t part_id;
        npart_in >> part_id >> std::ws;
        nodes_partition_ids.push_back(part_id);
    }
    npart_in.close();
    //TEST  std::remove(fname_parts.c_str());

    // -- Open files for output
    //
    // for ASCII output
    std::fstream os_subd;
    std::fstream os_subd_head;
    std::fstream os_subd_node;
    // for binary output
    FILE *of_bin_cfg = 0;   // File to output headers
    FILE *of_bin_nod = 0;   // File to output an array of all nodes
    FILE *of_bin_ele = 0;   // File to output an array of all non-ghost elements
    FILE *of_bin_ele_g = 0; // File to output an array of all ghost elements

    if (output_binary)
    {
        std::string fname = file_name_base + "_partitioned_msh_cfg" + npartitions_str + ".bin";
        of_bin_cfg = fopen(fname.c_str(), "wb");

        fname = file_name_base + "_partitioned_msh_ele" + npartitions_str + ".bin";
        of_bin_ele = fopen(fname.c_str(), "wb");

        fname = file_name_base + "_partitioned_msh_ele_g" + npartitions_str + ".bin";
        of_bin_ele_g = fopen(fname.c_str(), "wb");

        fname = file_name_base + "_partitioned_msh_nod" + npartitions_str + ".bin";
        of_bin_nod = fopen(fname.c_str(), "wb");
    }
    else
    {
        std::string fname = file_name_base + "_partitioned_cfg" + npartitions_str + ".msh";
        os_subd_head.open(fname.c_str(), std::ios::out|std::ios::trunc );
        const std::string mesh_info = "Subdomain mesh "
                                      "(Number of non-ghost nodes;  Number of non ghost base nodes; Number of elements; "
                                      "Number of ghost elements; Number of base nodes; Number of nodes) "
                                      "Number of the base nodes of the global mesh; Number of the nodes of the global mesh; "
                                      "Number of integers to define elements; Number of integers to define ghost elements; "
                                      "Reserved number.";
        os_subd_head << mesh_info << std::endl;
        os_subd_head << npartitions << std::endl;

        fname = file_name_base + "_partitioned_elems_" + npartitions_str + ".msh";
        os_subd.open(fname.c_str(), std::ios::out|std::ios::trunc );

        fname = file_name_base + "_partitioned_nodes_" + npartitions_str + ".msh";
        os_subd_node.open(fname.c_str(), std::ios::out|std::ios::trunc );

        std::setw(14);
        os_subd.precision(14);
        os_subd.setf(std::ios::scientific);
    }

    PetscInt global_node_id_offset = 0;
    std::vector<bool> nodes_part_reserved;
    nodes_part_reserved.reserve(nnodes);
    std::vector<PetscInt> node_global_ids;
    node_global_ids.reserve(nnodes);
    for (std::size_t i=0; i<nnodes; i++)
    {
        nodes_part_reserved.push_back(false);
        node_global_ids.push_back(0);
    }
    std::vector<bool> nodes_reseved(nnodes);
    std::vector<bool> elems_reseved(_elements.size());

    std::vector<MeshLib::Node*> partition_nodes;  // non-ghost nodes of a partition
    std::vector<std::size_t> partition_start_node_id(npartitions);
    std::vector<std::size_t> partition_end_node_id(npartitions);
    std::vector<std::size_t> end_non_ghost_node_id(npartitions);

    const PetscInt num_headers = 14;
    PetscInt head[num_headers];
    // node rank offset
    head[10] = 0;
    // element rank offset
    head[11] = 0;
    // ghost element rank offset
    head[12] = 0;
    // Reserved
    head[13] = 0;

    for (std::size_t ipart=0; ipart<npartitions; ipart++)
    {
        INFO("Processing partition: %d", ipart);

        partition_start_node_id[ipart] = partition_nodes.size();

        std::fill(nodes_reseved.begin(), nodes_reseved.end(), false);

        // Find non-ghost nodes in this partition
        for (std::size_t i=0; i< nnodes; i++)
        {
            if ( (nodes_partition_ids[i] == ipart) && (!nodes_part_reserved[i]) )
            {
                partition_nodes.push_back(_nodes[i]);
                nodes_part_reserved[i] = true;
                nodes_reseved[i] = true;
            }
        }
        end_non_ghost_node_id[ipart] = partition_nodes.size();

        // Find elements in this partition
        std::vector<const Element*> partition_regular_elements;
        std::vector<const Element*> partition_ghost_elements;
        std::fill(elems_reseved.begin(), elems_reseved.end(), false);
        // For ghost element
        std::vector< std::vector<unsigned> >
        elems_non_ghost_nodes_local_ids(_elements.size());

        const PetscInt num_non_ghost_nodes =  end_non_ghost_node_id[ipart]
                                           - partition_start_node_id[ipart];
        for (PetscInt i=0; i<num_non_ghost_nodes; i++)
        {
            const unsigned node_id = partition_nodes[ i +  partition_start_node_id[ipart] ]->getID();

            // Search the elements connected to this nodes
            std::vector<Element*> const& connected_elems = _nodes[node_id]->getElements();

            for (const auto elem : connected_elems)
            {
                // If checked
                if (elems_reseved[elem->getID()])
                    continue;

                std::vector<unsigned> nonghost_nodes_local_ids;
                unsigned ghost_counter = 0;
                for (unsigned kk=0; kk<elem->getNNodes(); kk++)
                {
                    if (nodes_reseved[elem->getNodeIndex(kk)])
                    {
                        nonghost_nodes_local_ids.push_back(kk);
                    }
                    else
                    {
                        ghost_counter++;
                    }
                }

                // All nodes of this element are inside this partition
                if (ghost_counter == 0)
                {
                    partition_regular_elements.push_back(elem);
                }
                else // ghost element
                {
                    partition_ghost_elements.push_back(elem);
                    elems_non_ghost_nodes_local_ids[elem->getID()]
                        = std::move(nonghost_nodes_local_ids);
                }

                elems_reseved[elem->getID()] = true;
            }
        }

        // Add ghost nodes in ghost elements to the node vector of this partition
        // nodes_reseved is reused without cleaning
        // Mark the non-ghost nodes for each ghost element
        for (const auto ghost_elem : partition_ghost_elements)
        {
            for (unsigned k=0; k<ghost_elem->getNNodes(); k++)
                nodes_reseved[ghost_elem->getNodeIndex(k)] = false;

            // Mark non-ghost nodes
            std::vector<unsigned>& local_node_ids = elems_non_ghost_nodes_local_ids[ghost_elem->getID()];
            for (std::size_t k=0; k<local_node_ids.size(); k++)
                nodes_reseved[ghost_elem->getNodeIndex(local_node_ids[k])] = true;
        }
        // Add the ghost nodes to the node vector of this partition
        for (const auto ghost_elem : partition_ghost_elements)
        {
            for (unsigned k=0; k<ghost_elem->getNNodes(); k++)
            {
                const unsigned node_id = ghost_elem->getNodeIndex(k);
                if (nodes_reseved[node_id])
                    continue;
                nodes_reseved[node_id] = true;
                partition_nodes.push_back(ghost_elem->getNode(k));
            }
        }

        partition_end_node_id[ipart] = partition_nodes.size();
        // Renumber
        PetscInt local_id = 0;

        std::vector<unsigned> nodes_local_ids_partition;
        nodes_local_ids_partition.reserve(nnodes);
        for (std::size_t i=0; i<nnodes; i++)
        {
            nodes_local_ids_partition.push_back(0);
        }

        for (std::size_t i = partition_start_node_id[ipart];
                i < end_non_ghost_node_id[ipart]; i++)
        {
            const unsigned node_id = partition_nodes[i]->getID();
            node_global_ids[node_id] = global_node_id_offset;
            nodes_local_ids_partition[node_id] = local_id;
            local_id++;
            global_node_id_offset++;
        }
        // Assign local IDs to ghost nodes
        for (std::size_t i = end_non_ghost_node_id[ipart];
                i < partition_end_node_id[ipart]; i++)
        {
            nodes_local_ids_partition[partition_nodes[i]->getID()] = local_id;
            local_id++;
        }

        // Count the number of integer variables used to define non-ghost elements
        const PetscInt num_non_ghost_elems = partition_regular_elements.size();
        PetscInt nmb_element_idxs = 3 * num_non_ghost_elems;
        for (PetscInt j=0; j<num_non_ghost_elems; j++)
        {
            nmb_element_idxs += partition_regular_elements[j]->getNNodes();
        }
        const PetscInt num_ghost_elems = partition_ghost_elements.size();
        PetscInt nmb_element_idxs_g = 3 * num_ghost_elems;
        for (PetscInt j=0; j<num_ghost_elems; j++)
        {
            nmb_element_idxs_g += partition_ghost_elements[j]->getNNodes();
        }

        std::string ipart_str = std::to_string(ipart);

        // Number of active elements
        const PetscInt offset_e = num_non_ghost_elems + nmb_element_idxs;
        const PetscInt offset_e_g = num_ghost_elems + nmb_element_idxs_g;

        // Reserved for nodes for high order interpolation
        const PetscInt num_non_ghost_nodes_extra = 0;

        const PetscInt num_nodes_linear_element =   partition_end_node_id[ipart]
                - partition_start_node_id[ipart];
        const PetscInt num_nodes_quadratic_element = 0;

        const PetscInt num_nodes_global_mesh = nnodes;
        // Reserved for nodes for high order interpolation
        const PetscInt num_nodes_global_mesh_extra = nnodes;

        // Write configuration data and elements of this partition
        head[0] = num_nodes_linear_element + num_nodes_quadratic_element;
        head[1] = num_nodes_linear_element;
        head[2] = num_non_ghost_elems;
        head[3] = num_ghost_elems;
        head[4] = num_non_ghost_nodes;
        head[5] = num_non_ghost_nodes + num_non_ghost_nodes_extra; // active nodes
        head[6] = num_nodes_global_mesh;
        head[7] = num_nodes_global_mesh_extra;
        head[8] = nmb_element_idxs;
        head[9] = nmb_element_idxs_g;
        if (output_binary)
        {
            fwrite(head, 1, num_headers * sizeof(PetscInt), of_bin_cfg);

            head[10] += head[0] * sizeof(NodeStruct);
            head[11] += offset_e * sizeof(PetscInt);
            head[12] += offset_e_g * sizeof(PetscInt);

            // Write non-ghost elements
            std::vector<PetscInt> ele_info(offset_e);

            PetscInt counter = num_non_ghost_elems;

            for (PetscInt j = 0; j < num_non_ghost_elems; j++)
            {
                ele_info[j] = counter;
                getElementIntegerVariables(*partition_regular_elements[j], nodes_local_ids_partition,
                                           ele_info, counter);
            }
            fwrite(&ele_info[0], 1, (offset_e) * sizeof(PetscInt), of_bin_ele);

            // Write ghost elements
            ele_info.resize(offset_e_g);
            counter = num_ghost_elems;
            for (PetscInt j = 0; j < num_ghost_elems; j++)
            {
                ele_info[j] = counter;
                getElementIntegerVariables(*partition_ghost_elements[j], nodes_local_ids_partition,
                                           ele_info, counter);
            }
            fwrite(&ele_info[0], 1, (offset_e_g) * sizeof(PetscInt), of_bin_ele_g);
        }
        else
        {

            for (int j=0; j<10; j++)
            {
                os_subd_head << head[j] << " ";
            }
            os_subd_head << " 0 " << std::endl;

            for (const auto elem : partition_regular_elements)
            {
                writeLocalElementNodeIndicies(os_subd, *elem, nodes_local_ids_partition);
            }

            for (const auto elem : partition_ghost_elements)
            {
                writeLocalElementNodeIndicies(os_subd, *elem, nodes_local_ids_partition);
            }
            os_subd << std::endl;
        }
    } // end of for(unsigned ipart=0; ipart<npartitions; ipart++) for partitioning of nodes, elements

    nodes_partition_ids.clear();
    nodes_part_reserved.clear();
    nodes_reseved.clear();
    elems_reseved.clear();

    // Write nodes
    if (output_binary)
    {
        fclose(of_bin_cfg);
        fclose(of_bin_ele);
        fclose(of_bin_ele_g);

        for (std::size_t ipart=0; ipart<npartitions; ipart++)
        {
            const size_t nnodes_part = partition_end_node_id[ipart] - partition_start_node_id[ipart];
            std::vector<NodeStruct> nodes_buffer;
            nodes_buffer.reserve(nnodes_part);

            for (std::size_t i=partition_start_node_id[ipart]; i<partition_end_node_id[ipart]; i++)
            {
                double const* coords = partition_nodes[i]->getCoords();
                NodeStruct node_struct;
                node_struct.id = node_global_ids[partition_nodes[i]->getID()];
                node_struct.x = coords[0];
                node_struct.y = coords[1];
                node_struct.z = coords[2];
                nodes_buffer.emplace_back(node_struct);
            }
            fwrite(&nodes_buffer[0], sizeof(NodeStruct), nnodes_part, of_bin_nod);
        }
        fclose(of_bin_nod);
    }
    else
    {
        os_subd_head.close();
        os_subd.close();
        std::setw(14);
        os_subd_node.precision(14);
        os_subd_node.setf(std::ios::scientific);

        for (std::size_t ipart=0; ipart<npartitions; ipart++)
        {
            for (std::size_t i=partition_start_node_id[ipart]; i<partition_end_node_id[ipart]; i++)
            {
                double const* coords = partition_nodes[i]->getCoords();
                os_subd_node << node_global_ids[partition_nodes[i]->getID()] << " "
                             << coords[0] << " " << coords[1] << " " << coords[2]<<"\n";
            }
            os_subd_node << std::endl;
        }
        os_subd_node.close();
    }

    // Update the node IDs
    for (std::size_t i=0; i<nnodes; i++)
    {
        _nodes[i]->setID(node_global_ids[i]);
    }
    // sort
    std::sort(_nodes.begin(), _nodes.end(),
              [](const MeshLib::Node* a, const MeshLib::Node* b)
                {
                   return a->getID() < b->getID();
                }
             );
}

ElementType MeshPartitioning::getElementType(const Element& elem)
{
    switch ( elem.getCellType() )
    {
        case MeshLib::CellType::LINE2:
            return ElementType::LINE2;
        case MeshLib::CellType::QUAD4:
            return ElementType::QUAD4;
        case MeshLib::CellType::HEX8:
            return ElementType::HEX8;
        case MeshLib::CellType::TRI3:
            return ElementType::TRI3;
        case MeshLib::CellType::PYRAMID5:
            return ElementType::PYRAMID5;
        case MeshLib::CellType::PRISM6:
            return ElementType::PRISM6;
        default:
            ERR("Invalid element type in element %d", elem.getID());
            abort();
    }
}

void MeshPartitioning::getElementIntegerVariables(const Element& elem,
        const std::vector<unsigned>& local_node_ids,
        std::vector<PetscInt>& elem_info,
        PetscInt& counter)
{
    unsigned mat_id = 0; // Materical ID to be set from the mesh data
    const PetscInt nn = elem.getNNodes();;
    elem_info[counter++] = mat_id;
    elem_info[counter++] = static_cast<unsigned>(getElementType(elem)) + 1;
    elem_info[counter++] = nn;

    for(PetscInt i=0; i<nn; i++)
    {
        elem_info[counter++] = local_node_ids[elem.getNodeIndex(i)];
    }
}

void MeshPartitioning::writeLocalElementNodeIndicies(std::ostream& os, const Element& elem,
        const std::vector<unsigned>& local_node_ids)
{
    unsigned mat_id = 0; // Materical ID to be set from the mesh data
    os << mat_id << " " << static_cast<unsigned>(getElementType(elem)) + 1 << " "
       << elem.getNNodes() <<" ";
    for(unsigned i=0; i<elem.getNNodes(); i++)
    {
        os << local_node_ids[elem.getNodeIndex(i)] << " ";
    }
    os << "\n";
}

}   // namespace MeshLib

