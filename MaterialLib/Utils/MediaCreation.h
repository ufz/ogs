// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

/// Parses a comma separated list of integers and/or ranges.
/// Such lists occur in the medium definition in the OGS prj file.
/// Range syntax is supported with colon separator, e.g., "1:5" expands to
/// "1,2,3,4,5". Multiple entries are separated by comma, e.g., "-1:2,5,7:9"
/// expands to "-1,0,1,2,5,7,8,9".
/// Error messages in this function refer to this specific purpose.
std::vector<int> splitMaterialIdString(std::string const& material_id_string);

/// Parses a comma separated list of integers or "*" string.
/// Such lists occur in the medium definition in the OGS prj file.
/// Range syntax is supported with colon separator, e.g., "1:5" expands to
/// "1,2,3,4,5". For the "*" input a vector of all (unique) material ids is
/// returned. Error messages in this function refer to this specific purpose.
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
