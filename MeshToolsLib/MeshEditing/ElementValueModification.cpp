/**
 * \file
 * \author Karsten Rink
 * \date   2013-04-04
 * \brief  Implementation of the ElementValueModification class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementValueModification.h"

#include <algorithm>
#include <range/v3/algorithm/fill.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>

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

    auto const old_values_filter = ranges::views::filter(
        [&old_value](auto const& v) { return v == old_value; });
    ranges::fill(*property_value_vec | old_values_filter, new_value);

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

    property_value_vector->assign(
        *property_value_vector |
        ranges::views::transform([&](auto const v)
                                 { return reverse_mapping[v]; }));

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

    auto const element_type_filter =
        ranges::views::filter([&](MeshLib::Element const* const e)
                              { return e->getGeomType() == ele_type; });

    auto selected_element_ids =
        mesh.getElements() | element_type_filter | MeshLib::views::ids;

    ranges::fill(
        selected_element_ids |
            ranges::views::transform([&](std::size_t const k) -> auto&
                                     { return (*property_value_vector)[k]; }),
        new_value);

    return ranges::distance(selected_element_ids);
}

}  // namespace MeshToolsLib
