/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Implementation of the MeshInformation class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshInformation.h"
#include "Elements/Element.h"

namespace MeshLib
{

const GeoLib::AABB MeshInformation::getBoundingBox(const MeshLib::Mesh &mesh)
{
    const std::vector<MeshLib::Node*> &nodes (mesh.getNodes());
    return GeoLib::AABB(nodes.begin(), nodes.end());
}

const std::array<unsigned, 7> MeshInformation::getNumberOfElementTypes(const MeshLib::Mesh &mesh)
{
    std::array<unsigned, 7> n_element_types = {{0, 0, 0, 0, 0, 0, 0}};
    const std::vector<MeshLib::Element*> &elements (mesh.getElements());
    for (auto it = elements.begin(); it != elements.end(); ++it)
    {
        MeshElemType t = (*it)->getGeomType();
        if (t == MeshElemType::LINE) n_element_types[0]++;
        if (t == MeshElemType::TRIANGLE) n_element_types[1]++;
        if (t == MeshElemType::QUAD) n_element_types[2]++;
        if (t == MeshElemType::TETRAHEDRON) n_element_types[3]++;
        if (t == MeshElemType::HEXAHEDRON) n_element_types[4]++;
        if (t == MeshElemType::PYRAMID) n_element_types[5]++;
        if (t == MeshElemType::PRISM) n_element_types[6]++;
    }
    return n_element_types;
}

} //end MeshLib
