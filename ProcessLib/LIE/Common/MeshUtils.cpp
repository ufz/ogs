/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"

#include "BaseLib/Algorithm.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
namespace
{
// A class to check whether a node is located on a crack tip with
// the following conditions:
// - the number of connected fracture elements is one
// - the node is not located on a domain boundary
class IsCrackTip
{
public:
    explicit IsCrackTip(MeshLib::Mesh const& mesh)
        : _mesh(mesh), _fracture_element_dim(mesh.getDimension() - 1)
    {
        _is_internal_node.resize(mesh.getNumberOfNodes(), true);

        MeshLib::NodeSearch nodeSearch(mesh);
        nodeSearch.searchBoundaryNodes();
        for (auto i : nodeSearch.getSearchedNodeIDs())
            _is_internal_node[i] = false;
    }

    bool operator()(MeshLib::Node const& node) const
    {
        if (!_is_internal_node[node.getID()] || !_mesh.isBaseNode(node.getID()))
            return false;

        unsigned n_connected_fracture_elements = 0;
        for (MeshLib::Element const* e : node.getElements())
            if (e->getDimension() == _fracture_element_dim)
                n_connected_fracture_elements++;
        assert(n_connected_fracture_elements > 0);

        return (n_connected_fracture_elements == 1);
    }

private:
    MeshLib::Mesh const& _mesh;
    unsigned const _fracture_element_dim;
    std::vector<bool> _is_internal_node;
};


void findFracutreIntersections(
    MeshLib::Mesh const& mesh,
    std::vector<int> const& vec_fracture_mat_IDs,
    std::vector<std::vector<MeshLib::Node*>> const& vec_fracture_nodes,
    std::vector<std::vector<MeshLib::Element*>>& intersected_fracture_elements,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_branch_nodeID_matIDs,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_junction_nodeID_matIDs
    )
{
    auto const n_fractures = vec_fracture_mat_IDs.size();
    std::map<unsigned, unsigned> matID_to_fid;
    for (unsigned i=0; i<n_fractures; i++)
        matID_to_fid[vec_fracture_mat_IDs[i]] = i;

    // make a vector all fracture nodes
    std::vector<std::size_t> all_fracture_nodes;
    for (auto& vec : vec_fracture_nodes)
        for (auto* node : vec)
            all_fracture_nodes.push_back(node->getID());

    // create a table of a node id and connected material IDs
    std::map<std::size_t, std::vector<std::size_t>> frac_nodeID_to_matIDs;
    for (unsigned i=0; i<n_fractures; i++)
        for (auto* node : vec_fracture_nodes[i])
            frac_nodeID_to_matIDs[node->getID()].push_back(vec_fracture_mat_IDs[i]);

    auto opt_material_ids(
        mesh.getProperties().getPropertyVector<int>("MaterialIDs"));

    // find branch/junction nodes which connect to multiple fractures
    intersected_fracture_elements.resize(n_fractures);
    for (auto entry : frac_nodeID_to_matIDs)
    {
        auto nodeID = entry.first;
        auto const* node = mesh.getNode(entry.first);
        auto const& matIDs = entry.second;
        if (matIDs.size() < 2)
            continue; // no intersection

        std::vector<MeshLib::Element*> conn_fracture_elements;
        {
            for (auto const* e : node->getElements())
                if (e->getDimension() == (mesh.getDimension()-1))
                    conn_fracture_elements.push_back(
                        const_cast<MeshLib::Element*>(e));
        }

        std::map<int,int> vec_matID_counts;
        {
            for (auto matid : matIDs)
                vec_matID_counts[matid] = 0;

            for (auto const* e : conn_fracture_elements)
            {
                auto matid = (*opt_material_ids)[e->getID()];
                vec_matID_counts[matid]++;
            }
        }

        for (auto matid : matIDs)
        {
            auto fid =  matID_to_fid[matid];
            for (auto* e : conn_fracture_elements)
            {
                auto e_matid = (*opt_material_ids)[e->getID()];
                if (matID_to_fid[e_matid] != fid) // only slave elements
                    intersected_fracture_elements[fid].push_back(e);
            }
        }

        bool isBranch = false;
        {
            for (auto entry : vec_matID_counts)
            {
                auto count = entry.second;
                if (count%2==1) {
                    isBranch = true;
                    break;
                }
            }
        }

        if (isBranch)
        {
            std::vector<int> matIDs(2);
            for (auto entry : vec_matID_counts)
            {
                auto matid = entry.first;
                auto count = entry.second;
                if (count%2==0) {
                    matIDs[0] = matid; // master
                } else {
                    matIDs[1] = matid; // slave
                }
            }
            vec_branch_nodeID_matIDs.push_back(std::make_pair(nodeID,matIDs));

        } else {
            std::vector<int> matIDs(2);
            matIDs[0] = std::min(vec_matID_counts.begin()->first, vec_matID_counts.rbegin()->first);
            matIDs[1] = std::max(vec_matID_counts.begin()->first, vec_matID_counts.rbegin()->first);
            vec_junction_nodeID_matIDs.push_back(std::make_pair(nodeID,matIDs));
        }
    }

    for (auto& eles : intersected_fracture_elements)
        BaseLib::makeVectorUnique(eles);

     DBUG("-> found %d branchs and %d junctions", vec_branch_nodeID_matIDs.size(), vec_junction_nodeID_matIDs.size());
}

}  // namespace

void getFractureMatrixDataInMesh(
    MeshLib::Mesh const& mesh,
    std::vector<MeshLib::Element*>& vec_matrix_elements,
    std::vector<int>& vec_fracture_mat_IDs,
    std::vector<std::vector<MeshLib::Element*>>& vec_fracture_elements,
    std::vector<std::vector<MeshLib::Element*>>& vec_fracture_matrix_elements,
    std::vector<std::vector<MeshLib::Node*>>& vec_fracture_nodes,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_branch_nodeID_matIDs,
    std::vector<std::pair<std::size_t,std::vector<int>>>& vec_junction_nodeID_matIDs)
{
    IsCrackTip isCrackTip(mesh);

    // get vectors of matrix elements and fracture elements
    vec_matrix_elements.reserve(mesh.getNumberOfElements());
    std::vector<MeshLib::Element*> all_fracture_elements;
    for (MeshLib::Element* e : mesh.getElements())
    {
        if (e->getDimension() == mesh.getDimension())
            vec_matrix_elements.push_back(e);
        else
            all_fracture_elements.push_back(e);
    }
    DBUG("-> found total %d matrix elements and %d fracture elements",
         vec_matrix_elements.size(), all_fracture_elements.size());

    // get fracture material IDs
    auto const material_ids = materialIDs(mesh);
    if (!material_ids)
    {
        OGS_FATAL("Could not access MaterialIDs property from mesh.");
    }
    for (MeshLib::Element* e : all_fracture_elements)
        vec_fracture_mat_IDs.push_back((*material_ids)[e->getID()]);
    BaseLib::makeVectorUnique(vec_fracture_mat_IDs);
    DBUG("-> found %d fracture material groups", vec_fracture_mat_IDs.size());

    // create a vector of fracture elements for each material
    vec_fracture_elements.resize(vec_fracture_mat_IDs.size());
    for (unsigned frac_id = 0; frac_id < vec_fracture_mat_IDs.size(); frac_id++)
    {
        const auto frac_mat_id = vec_fracture_mat_IDs[frac_id];
        std::vector<MeshLib::Element*>& vec_elements =
            vec_fracture_elements[frac_id];
        std::copy_if(all_fracture_elements.begin(), all_fracture_elements.end(),
                     std::back_inserter(vec_elements),
                     [&](MeshLib::Element* e) {
                         return (*material_ids)[e->getID()] == frac_mat_id;
                     });
        DBUG("-> found %d elements on the fracture %d", vec_elements.size(),
             frac_id);
    }

    // get a vector of fracture nodes for each material
    vec_fracture_nodes.resize(vec_fracture_mat_IDs.size());
    for (unsigned frac_id = 0; frac_id < vec_fracture_mat_IDs.size(); frac_id++)
    {
        std::vector<MeshLib::Node*>& vec_nodes = vec_fracture_nodes[frac_id];
        for (MeshLib::Element* e : vec_fracture_elements[frac_id])
        {
            for (unsigned i = 0; i < e->getNumberOfNodes(); i++)
            {
                if (isCrackTip(*e->getNode(i)))
                    continue;
                vec_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
            }
        }
        BaseLib::makeVectorUnique(
            vec_nodes, [](MeshLib::Node* node1, MeshLib::Node* node2) {
                return node1->getID() < node2->getID();
            });
        DBUG("-> found %d nodes on the fracture %d", vec_nodes.size(), frac_id);
    }

    // find branch/junction nodes which connect to multiple fractures
    std::vector<std::vector<MeshLib::Element*>> intersected_fracture_elements;
    findFracutreIntersections(
        mesh, vec_fracture_mat_IDs, vec_fracture_nodes,
        intersected_fracture_elements,
        vec_branch_nodeID_matIDs, vec_junction_nodeID_matIDs);

    // create a vector fracture elements and connected matrix elements,
    // which are passed to a DoF table
    for (unsigned fid = 0; fid<vec_fracture_elements.size(); fid++)
    {
        auto const& fracture_elements = vec_fracture_elements[fid];
        std::vector<MeshLib::Element*> vec_ele;
        // first, collect matrix elements
        for (MeshLib::Element* e : fracture_elements)
        {
            // it is sufficient to iterate over base nodes, because they are
            // already connected to all neighbours
            for (unsigned i = 0; i < e->getNumberOfBaseNodes(); i++)
            {
                MeshLib::Node const* node = e->getNode(i);
                if (isCrackTip(*node))
                    continue;
                for (unsigned j = 0; j < node->getNumberOfElements(); j++)
                {
                    // only matrix elements
                    if (node->getElement(j)->getDimension() <
                        mesh.getDimension())
                        continue;
                    vec_ele.push_back(
                        const_cast<MeshLib::Element*>(node->getElement(j)));
                }
            }
        }
        BaseLib::makeVectorUnique(
            vec_ele, [](MeshLib::Element* e1, MeshLib::Element* e2) {
                return e1->getID() < e2->getID();
            });

        // second, append fracture elements
        std::copy(fracture_elements.begin(), fracture_elements.end(),
                  std::back_inserter(vec_ele));
        // thirdly, append intersected fracture elements
        std::copy(intersected_fracture_elements[fid].begin(), intersected_fracture_elements[fid].end(),
                  std::back_inserter(vec_ele));

        vec_fracture_matrix_elements.push_back(vec_ele);
    }
}

}  // namespace LIE
}  // namespace ProcessLib
