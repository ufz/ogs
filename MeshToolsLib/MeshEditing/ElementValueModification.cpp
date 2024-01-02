/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the ElementValueModification class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementValueModification.h"

#include <algorithm>

#include "BaseLib/Logging.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

namespace MeshToolsLib
{
bool ElementValueModification::replace(MeshLib::Mesh& mesh,
                                       std::string const& property_name,
                                       int const old_value, int const new_value,
                                       bool replace_if_exists)
{
    MeshLib::PropertyVector<int>* property_value_vec = nullptr;
    try
    {
        property_value_vec = mesh.getProperties().getPropertyVector<int>(
            property_name, MeshLib::MeshItemType::Cell, 1);
    }
    catch (std::runtime_error const& e)
    {
        ERR("{:s}", e.what());
        return false;
    }

    const std::size_t n_property_tuples(
        property_value_vec->getNumberOfTuples());

    if (!replace_if_exists)
    {
        for (std::size_t i = 0; i < n_property_tuples; ++i)
        {
            if ((*property_value_vec)[i] == new_value)
            {
                WARN(
                    "ElementValueModification::replaceElementValue() - "
                    "Replacement value '{:d}' is already taken, no changes "
                    "have been made.",
                    new_value);
                return false;
            }
        }
    }

    for (std::size_t i = 0; i < n_property_tuples; ++i)
    {
        if ((*property_value_vec)[i] == old_value)
        {
            (*property_value_vec)[i] = new_value;
        }
    }

    return true;
}

bool ElementValueModification::replace(MeshLib::Mesh& mesh, int const old_value,
                                       int const new_value,
                                       bool replace_if_exists)
{
    return replace(mesh, "MaterialIDs", old_value, new_value,
                   replace_if_exists);
}

std::size_t ElementValueModification::condense(MeshLib::Mesh& mesh)
{
    MeshLib::PropertyVector<int>* property_value_vector = nullptr;
    try
    {
        property_value_vector = mesh.getProperties().getPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    }
    catch (std::runtime_error const& e)
    {
        ERR("{:s}", e.what());
        return 0;
    }

    std::vector<int> value_mapping(
        getSortedPropertyValues(*property_value_vector));

    std::vector<int> reverse_mapping(value_mapping.back() + 1, 0);
    std::size_t const nValues(value_mapping.size());
    for (std::size_t i = 0; i < nValues; ++i)
    {
        reverse_mapping[value_mapping[i]] = i;
    }

    std::size_t const n_property_values(property_value_vector->size());
    for (std::size_t i = 0; i < n_property_values; ++i)
    {
        (*property_value_vector)[i] =
            reverse_mapping[(*property_value_vector)[i]];
    }

    return nValues;
}

std::size_t ElementValueModification::setByElementType(
    MeshLib::Mesh& mesh, MeshLib::MeshElemType ele_type, int const new_value)
{
    MeshLib::PropertyVector<int>* property_value_vector = nullptr;
    try
    {
        property_value_vector = mesh.getProperties().getPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    }
    catch (std::runtime_error const& e)
    {
        ERR("{:s}", e.what());
        return 0;
    }

    std::vector<MeshLib::Element*> const& elements(mesh.getElements());
    std::size_t cnt(0);
    for (std::size_t k(0); k < elements.size(); k++)
    {
        if (elements[k]->getGeomType() != ele_type)
        {
            continue;
        }
        (*property_value_vector)[k] = new_value;
        cnt++;
    }

    return cnt;
}

}  // namespace MeshToolsLib
