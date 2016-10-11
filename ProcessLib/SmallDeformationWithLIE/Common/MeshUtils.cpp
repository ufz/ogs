/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshUtils.h"


namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

void getFractureMatrixDataInMesh(
        MeshLib::Mesh const& mesh,
        std::vector<MeshLib::Element*>& vec_matrix_elements,
        std::vector<MeshLib::Element*>& vec_fracture_elements,
        std::vector<MeshLib::Element*>& vec_fracture_matrix_elements,
        std::vector<MeshLib::Node*>& vec_fracture_nodes
        )
{
    // get vectors of matrix elements and fracture elements
    vec_matrix_elements.reserve(mesh.getNumberOfElements());
    for (MeshLib::Element* e : mesh.getElements())
    {
        if (e->getDimension() == mesh.getDimension())
            vec_matrix_elements.push_back(e);
        else
            vec_fracture_elements.push_back(e);
    }
    DBUG("-> found total %d matrix elements and %d fracture elements",
         vec_matrix_elements.size(), vec_fracture_elements.size());

    // get a vector of fracture nodes
    for (MeshLib::Element* e : vec_fracture_elements)
    {
        for (unsigned i=0; i<e->getNumberOfNodes(); i++)
        {
            vec_fracture_nodes.push_back(const_cast<MeshLib::Node*>(e->getNode(i)));
        }
    }
    std::sort(vec_fracture_nodes.begin(), vec_fracture_nodes.end(),
        [](MeshLib::Node* node1, MeshLib::Node* node2) { return (node1->getID() < node2->getID()); }
        );
    vec_fracture_nodes.erase(
                std::unique(vec_fracture_nodes.begin(), vec_fracture_nodes.end()),
                vec_fracture_nodes.end());
    DBUG("-> found %d nodes on the fracture", vec_fracture_nodes.size());

    // create a vector fracture elements and connected matrix elements,
    // which are passed to a DoF table
    // first, collect matrix elements
    for (MeshLib::Element *e : vec_fracture_elements)
    {
        for (unsigned i=0; i<e->getNumberOfBaseNodes(); i++)
        {
            MeshLib::Node const* node = e->getNode(i);
            for (unsigned j=0; j<node->getNumberOfElements(); j++)
            {
                // only matrix elements
                if (node->getElement(j)->getDimension() == mesh.getDimension()-1)
                    continue;
                vec_fracture_matrix_elements.push_back(const_cast<MeshLib::Element*>(node->getElement(j)));
            }
        }
    }
    std::sort(vec_fracture_matrix_elements.begin(), vec_fracture_matrix_elements.end(),
        [](MeshLib::Element* p1, MeshLib::Element* p2) { return (p1->getID() < p2->getID()); }
        );
    vec_fracture_matrix_elements.erase(
                std::unique(vec_fracture_matrix_elements.begin(), vec_fracture_matrix_elements.end()),
                vec_fracture_matrix_elements.end());

    // second, append fracture elements
    vec_fracture_matrix_elements.insert(
                vec_fracture_matrix_elements.end(),
                vec_fracture_elements.begin(), vec_fracture_elements.end());
}

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib
