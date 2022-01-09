/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Implementation of the MeshInformation class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshInformation.h"

#include "Elements/Element.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

namespace
{
template <typename T>
void printBounds(MeshLib::PropertyVector<T> const& property)
{
    auto const bounds = MeshLib::MeshInformation::getValueBounds(property);
    if (!bounds.has_value())
    {
        INFO("\t{:s}: Could not get value bounds for property vector.",
             property.getPropertyName());
        return;
    }
    INFO("\t{:s}: ({:d} values) [{}, {}]", property.getPropertyName(),
         property.size(), bounds->first, bounds->second);
}
}  // namespace

namespace MeshLib
{
GeoLib::AABB MeshInformation::getBoundingBox(const MeshLib::Mesh& mesh)
{
    const std::vector<MeshLib::Node*>& nodes(mesh.getNodes());
    return GeoLib::AABB(nodes.begin(), nodes.end());
}

std::map<MeshElemType, unsigned> MeshInformation::getNumberOfElementTypes(
    const MeshLib::Mesh& mesh)
{
    std::map<MeshElemType, unsigned> n_element_types;
    const std::vector<MeshLib::Element*>& elements(mesh.getElements());
    for (auto element : elements)
    {
        MeshElemType t = element->getGeomType();
        n_element_types[t]++;
    }
    return n_element_types;
}

void MeshInformation::writeAllNumbersOfElementTypes(const MeshLib::Mesh& mesh)
{
    auto const& nr_ele_types =
        MeshLib::MeshInformation::getNumberOfElementTypes(mesh);

    INFO("Number of elements in the mesh:");
    for (auto entry : nr_ele_types)
    {
        INFO("\t{:s}s: {:d}",
             MeshLib::MeshElemType2String(
                 static_cast<MeshLib::MeshElemType>(entry.first)),
             entry.second);
    }
}

void MeshInformation::writePropertyVectorInformation(const MeshLib::Mesh& mesh)
{
    auto const& properties = mesh.getProperties();
    INFO("There are {:d} properties in the mesh:", properties.size());

    for (auto [name, property] : properties)
    {
        if (auto p = dynamic_cast<PropertyVector<double>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<float>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<int>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<long>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<long long>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p =
                     dynamic_cast<PropertyVector<unsigned long>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<unsigned long long>*>(
                     property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<std::size_t>*>(property))
        {
            printBounds(*p);
        }
        else if (auto p = dynamic_cast<PropertyVector<char>*>(property))
        {
            printBounds(*p);
        }
        else
        {
            INFO(
                "\t{:s}: Could not get value bounds for property vector of "
                "type '{:s}'.",
                name, typeid(*p).name());
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
