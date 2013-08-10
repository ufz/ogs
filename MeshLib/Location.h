/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>

#ifndef LOCATION_H_
#define LOCATION_H_

namespace MeshLib
{

enum class MeshItemType { Node, Edge, Face, Cell };

std::ostream& operator<<(std::ostream& os, MeshItemType const& t);

/// Spatial location description.
///
/// The spatial location is given by a mesh by its \c mesh_id, item's type (face,
/// cell, etc. see MeshItemType), and item's number by its \c item_id.
struct Location
{
    std::size_t          mesh_id;
    MeshItemType         item_type;
    std::size_t          item_id;

    Location(std::size_t meshid, MeshItemType itemtype, std::size_t itemid)
    : mesh_id(meshid), item_type(itemtype), item_id(itemid){};
};

/// Lexicographic order of Location.
bool operator<(const Location& left, const Location& right);

std::ostream& operator<<(std::ostream& os, Location const& l);

}   // namespace MeshLib

#endif  // LOCATION_H_
