/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Implementation of the MeshInformation class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshInformation.h"

#include "Elements/Element.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

namespace MeshLib
{
GeoLib::AABB MeshInformation::getBoundingBox(const MeshLib::Mesh& mesh)
{
    const std::vector<MeshLib::Node*>& nodes(mesh.getNodes());
    return GeoLib::AABB(nodes.begin(), nodes.end());
}

std::array<unsigned, 7> MeshInformation::getNumberOfElementTypes(
    const MeshLib::Mesh& mesh)
{
    std::array<unsigned, 7> n_element_types = {{0, 0, 0, 0, 0, 0, 0}};
    const std::vector<MeshLib::Element*>& elements(mesh.getElements());
    for (auto element : elements)
    {
        MeshElemType t = element->getGeomType();
        if (t == MeshElemType::LINE)
        {
            n_element_types[0]++;
        }
        if (t == MeshElemType::TRIANGLE)
        {
            n_element_types[1]++;
        }
        if (t == MeshElemType::QUAD)
        {
            n_element_types[2]++;
        }
        if (t == MeshElemType::TETRAHEDRON)
        {
            n_element_types[3]++;
        }
        if (t == MeshElemType::HEXAHEDRON)
        {
            n_element_types[4]++;
        }
        if (t == MeshElemType::PYRAMID)
        {
            n_element_types[5]++;
        }
        if (t == MeshElemType::PRISM)
        {
            n_element_types[6]++;
        }
    }
    return n_element_types;
}

void MeshInformation::writeAllNumbersOfElementTypes(const MeshLib::Mesh& mesh)
{
    auto const& nr_ele_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(mesh);
    const std::array<std::string, 7> element_names{
        "lines",     "triangles", "quads", "tetrahedra",
        "hexahedra", "pyramids",  "prisms"};
    for (unsigned int element_type = 0; element_type < nr_ele_types.size();
         element_type++)
    {
        INFO("\t{:d} {:s} ", nr_ele_types[element_type],
             element_names[element_type].c_str());
    }
}

void MeshInformation::writePropertyVectorInformation(const MeshLib::Mesh& mesh)
{
    std::vector<std::string> const& vec_names(
        mesh.getProperties().getPropertyVectorNames());
    INFO("There are {:d} properties in the mesh:", vec_names.size());
    for (const auto& vec_name : vec_names)
    {
        if (auto const vec_bounds =
                MeshLib::MeshInformation::getValueBounds<int>(mesh, vec_name))
        {
            INFO("\t{:s}: [{:d}, {:d}]", vec_name.c_str(), vec_bounds->first,
                 vec_bounds->second);
        }
        else if (auto const vec_bounds =
                     MeshLib::MeshInformation::getValueBounds<double>(mesh,
                                                                      vec_name))
        {
            INFO("\t{:s}: [{:g}, {:g}]", vec_name.c_str(), vec_bounds->first,
                 vec_bounds->second);
        }
        else
        {
            INFO("\t{:s}: Could not get value bounds for property vector.",
                 vec_name.c_str());
        }
    }
}

void MeshInformation::writeMeshValidationResults(MeshLib::Mesh& mesh)
{
    MeshLib::MeshValidation validation(mesh);

    unsigned const n_holes(MeshLib::MeshValidation::detectHoles(mesh));
    if (n_holes > 0)
    {
        INFO("{:d} hole(s) detected within the mesh", n_holes);
    }
    else
    {
        INFO("No holes found within the mesh.");
    }
}

}  // namespace MeshLib
