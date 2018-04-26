/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>
#include <vector>

#include "GeoLib/AABB.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib {

// forward declarations
class Mesh;
class Element;

/// Element search class
class ElementSearch final
{
public:
    explicit ElementSearch(const MeshLib::Mesh &mesh);

    /// return marked elements
    const std::vector<std::size_t>& getSearchedElementIDs() const { return _marked_elements; }

    /// @tparam PROPERTY_TYPE integral type of the property
    /// Different properties can be assigned to the elements of the mesh. These
    /// properties can be accessed by the name of the property. The method marks
    /// all elements of the mesh for the property \c property_name with  a
    /// property value equal to \c property_value.
    /// @param property_name the name of the property the searching/marking is
    /// based on
    /// @param property_value value required for the element to be marked
    /// @return The number of marked elements will be returned. The concrete
    /// element ids can be requested by getSearchedElementIDs().
    template <typename PROPERTY_TYPE>
    std::size_t searchByPropertyValue(std::string const& property_name,
                                      PROPERTY_TYPE const property_value)
    {
        return searchByPropertyValueRange<PROPERTY_TYPE>(
            property_name, property_value, property_value, false);
    }

    /// @tparam PROPERTY_TYPE integral type of the property
    /// Different properties can be assigned to the elements of the mesh. These
    /// properties can be accessed by the name of the property. The method marks
    /// all elements of the mesh for the property \c property_name with  a
    /// property value outside of the interval [min_property_value,
    /// max_property_value].
    /// @param property_name the name of the property the searching/marking is
    /// based on
    /// @param min_property_value minimum value of the given property for the
    /// element not to be marked
    /// @param max_property_value maximum value of the given property for the
    /// element not to be marked
    /// @param outside_of if true, all values outside of the given range or
    /// markedc, if false, all values inside the given range are marked
    /// @return The number of marked elements will be returned. The concrete
    /// element ids can be requested by getSearchedElementIDs().
    template <typename PROPERTY_TYPE>
    std::size_t searchByPropertyValueRange(
        std::string const& property_name,
        PROPERTY_TYPE const min_property_value,
        PROPERTY_TYPE const max_property_value,
        bool outside_of = true)
    {
        if (!_mesh.getProperties().existsPropertyVector<PROPERTY_TYPE>(property_name))
        {
            WARN("Property \"%s\" not found in mesh.", property_name.c_str());
            return 0;
        }
        auto const* const pv =
            _mesh.getProperties().getPropertyVector<PROPERTY_TYPE>(property_name);

        if (pv->getMeshItemType() != MeshLib::MeshItemType::Cell)
        {
            WARN("The property \"%s\" is not assigned to mesh elements.",
                 property_name.c_str());
            return 0;
        }

        std::vector<std::size_t> matchedIDs;

        if (outside_of)
        {
            for (std::size_t i(0); i < pv->getNumberOfTuples(); ++i)
            {
                if ((*pv)[i] < min_property_value || (*pv)[i] > max_property_value)
                    matchedIDs.push_back(i);
            }
        }
        else
        {
            for (std::size_t i(0); i < pv->getNumberOfTuples(); ++i)
            {
                if ((*pv)[i] >= min_property_value && (*pv)[i] <= max_property_value)
                    matchedIDs.push_back(i);
            }
        }
        updateUnion(matchedIDs);
        return matchedIDs.size();
    }

    /// Marks all elements of the given element type.
    std::size_t searchByElementType(MeshElemType eleType);

    /// Marks all elements with a volume smaller than eps.
    std::size_t searchByContent(double eps = std::numeric_limits<double>::epsilon());

    /// Marks all elements with at least one node outside the bounding box spanned by x1 and x2;
    std::size_t searchByBoundingBox(GeoLib::AABB const& aabb);

    /// Marks all elements connecting to any of the given nodes
    std::size_t searchByNodeIDs(const std::vector<std::size_t> &node_ids);

private:
    /// Updates the vector of marked elements with values from vec.
    void updateUnion(const std::vector<std::size_t> &vec);

    /// The mesh from which elements should be removed.
    const MeshLib::Mesh &_mesh;
    /// The vector of element indices that should be removed.
    std::vector<std::size_t> _marked_elements;
};

} // end namespace MeshLib
