/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENTSEARCH_H_
#define ELEMENTSEARCH_H_

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
    /// all elements of the mesh for the property \c property_name with the
    /// given property value \c property_value.
    /// @param property_value the value of the property the elements have to
    /// have to be marked
    /// @param property_name the name of the property the searching/marking is
    /// based on
    /// @return The number of marked elements will be returned. The concrete
    /// element ids can be requested by getSearchedElementIDs().
    template <typename PROPERTY_TYPE>
    std::size_t searchByPropertyValue(
        PROPERTY_TYPE const property_value,
        std::string const& property_name = "MaterialIDs")
    {
        auto const* const pv =
            _mesh.getProperties().getPropertyVector<PROPERTY_TYPE>(
                property_name);
        if (!pv)
        {
            WARN("Property \"%s\" not found in mesh.", property_name.c_str());
            return 0;
        }

        if (pv->getMeshItemType() != MeshLib::MeshItemType::Cell)
        {
            WARN("The property \"%s\" is not assigned to mesh elements.",
                 property_name.c_str());
            return 0;
        }

        std::vector<std::size_t> matchedIDs;
        for (std::size_t i(0); i < pv->getNumberOfTuples(); ++i) {
            if ((*pv)[i] == property_value)
                matchedIDs.push_back(i);
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

#endif //ELEMENTEXTRACTION_H
