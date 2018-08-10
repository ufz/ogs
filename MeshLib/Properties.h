/**
 * \file
 * \brief  Definition of the class Properties that implements a container of
 *         properties.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdlib>
#include <string>
#include <map>

#include <logog/include/logog.hpp>

#include "Location.h"

#include "PropertyVector.h"

namespace MeshLib
{

/// @brief Property manager on mesh items.
/// Class Properties manages scalar, vector or matrix properties. For instance
/// in groundwater flow porosity is a scalar property and permeabilty can be
/// stored as a tensor property. Properties are assigned to mesh items, i.e.
/// Node or Element objects. The createNewPropertyVector() method first creates a
/// PropertyVector of template type T (scalar, vector or matrix).
/// This class stores the PropertyVector, accessible by a combination of the
/// name and the type of the mesh item (Node or Element).
class Properties
{
public:
    /// Method creates a PropertyVector if a PropertyVector with the same name
    /// and the same type T was not already created before. In case there exists
    /// already such a PropertyVector a nullptr is returned.
    /// There are two versions of this method. This method is used when every
    /// mesh item at hand has its own property value, i.e. \f$n\f$ mesh item and
    /// \f$n\f$ different property values.
    /// The user has to ensure the correct usage of the vector later on.
    /// @tparam T type of the property value
    /// @param name the name of the property
    /// @param mesh_item_type for instance node or element assigned properties
    /// @param n_components number of components for each tuple
    /// @return A pointer to a PropertyVector on success and a nullptr
    /// otherwise.
    template <typename T>
    PropertyVector<T>* createNewPropertyVector(std::string const& name,
                                               MeshItemType mesh_item_type,
                                               std::size_t n_components = 1);

    /// Method creates a PropertyVector if a PropertyVector with the same name
    /// and the same type T was not already created before. In case there exists
    /// already such a PropertyVector a nullptr is returned.
    /// This method is used if only a small number of distinct property values
    /// in a property exist (e.g. mapping property groups to elements).
    /// In this case a mapping between mesh items and properties (stored
    /// on the heap), see the parameter item2group_mapping, is required.
    /// @tparam T type of the property value
    /// @param name the name of the property
    /// @param n_prop_groups number of distinct property groups
    /// @param item2group_mapping the mapping between mesh item and the property
    /// group
    /// @param mesh_item_type for instance node or element assigned properties
    /// @param n_components number of components for each tuple
    /// @return A pointer to a PropertyVector on success and a nullptr
    /// otherwise.
    template <typename T>
    PropertyVector<T>* createNewPropertyVector(
        std::string const& name,
        std::size_t n_prop_groups,
        std::vector<std::size_t> const& item2group_mapping,
        MeshItemType mesh_item_type,
        std::size_t n_components = 1);

    /// Checks if a property vector with given \c name and the given type
    /// exists.
    /// @param name name of the requested property vector
    template <typename T>
    bool existsPropertyVector(std::string const& name) const;

    /// Returns a property vector with given \c name or nullptr if no such
    /// property vector exists.
    template <typename T>
    PropertyVector<T> const* getPropertyVector(std::string const& name) const;

    /// Returns a property vector with given \c name or nullptr if no such
    /// property vector exists.
    template <typename T>
    PropertyVector<T>* getPropertyVector(std::string const& name);

    void removePropertyVector(std::string const& name);

    /// Check if a PropertyVector accessible by the name is already
    /// stored within the Properties object.
    /// @param name the name of the property (for instance porosity)
    bool hasPropertyVector(std::string const& name) const;

    std::vector<std::string> getPropertyVectorNames() const;
    std::vector<std::string> getPropertyVectorNames(
        MeshLib::MeshItemType t) const;

    /** copy all PropertyVector objects stored in the (internal) map but only
     * those nodes/elements of a PropertyVector whose ids are not in the vectors
     * exclude_*_ids.
     */
    Properties excludeCopyProperties(
        std::vector<std::size_t> const& exclude_elem_ids,
        std::vector<std::size_t> const& exclude_node_ids) const;

    /** copy all PropertyVector objects stored in the (internal) map but
     * PropertyVector objects with the given MeshItemType are excluded.
     */
    Properties excludeCopyProperties(
        std::vector<MeshItemType> const& exclude_mesh_item_types) const;

    Properties() = default;

    Properties(Properties const& properties);
    Properties(Properties&& properties) = default;
    Properties& operator=(Properties const& properties);
    Properties& operator=(Properties&& properties) = default;

    ~Properties();

private:
    /// A mapping from property's name to the stored object of any type.
    /// See addProperty() and getProperty() documentation.
    std::map<std::string, PropertyVectorBase*> _properties;
}; // end class

#include "Properties-impl.h"

} // end namespace MeshLib
