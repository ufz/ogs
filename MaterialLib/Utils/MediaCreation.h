/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/Error.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}  // namespace MeshLib

namespace MaterialLib
{

/// Parses a comma separated list of integers.
/// Such lists occur in the medium definition in the OGS prj file.
/// Error messages in this function refer to this specific purpose.
std::vector<int> splitMaterialIdString(std::string const& material_id_string);

/// Parses a comma separated list of integers or "*" string.
/// Such lists occur in the medium definition in the OGS prj file.
/// For the "*" input a vector of all (unique) material ids is returned.
/// Error messages in this function refer to this specific purpose.
std::vector<int> parseMaterialIdString(
    std::string const& material_id_string,
    MeshLib::PropertyVector<int> const* const material_ids);

/// Creates a new entry for the material id in the media map by either calling
/// the create_medium function and creating a new shared pointer, or by reusing
/// the existing shared pointer.
template <typename T, typename CreateMedium>
    requires std::convertible_to<
        decltype(std::declval<CreateMedium>()(std::declval<int>())),
        std::shared_ptr<T>>
void createMediumForId(int const id,
                       std::map<int, std::shared_ptr<T>>& media,
                       std::vector<int> const& material_ids_of_this_medium,
                       CreateMedium&& create_medium)
{
    if (media.find(id) != end(media))
    {
        OGS_FATAL(
            "Multiple media were specified for the same material id '{:d}'. "
            "Keep in mind, that if no material id is specified, it is assumed "
            "to be 0 by default.",
            id);
    }

    if (id == material_ids_of_this_medium[0])
    {
        media[id] = create_medium(id);
    }
    else
    {
        media[id] = media[material_ids_of_this_medium[0]];
    }
}

}  // namespace MaterialLib
