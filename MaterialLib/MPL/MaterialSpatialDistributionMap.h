// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <range/v3/view.hpp>
#include <vector>

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}  // namespace MeshLib

namespace MaterialPropertyLib
{
class Medium;

class MaterialSpatialDistributionMap
{
public:
    MaterialSpatialDistributionMap(
        std::map<int, std::shared_ptr<Medium>> const& media,
        MeshLib::PropertyVector<int> const* const material_ids)
        : media_(media), material_ids_(material_ids)
    {
    }

    auto media() const { return media_ | ranges::views::values; }

    Medium* getMedium(std::size_t element_id);
    Medium const* getMedium(std::size_t element_id) const;
    void checkElementHasMedium(std::size_t const element_id) const;

private:
    std::map<int, std::shared_ptr<Medium>> const& media_;
    MeshLib::PropertyVector<int> const* const material_ids_;
};
}  // namespace MaterialPropertyLib
