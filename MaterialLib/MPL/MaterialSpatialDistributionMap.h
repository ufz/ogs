/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <map>
#include <memory>
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

    Medium* getMedium(std::size_t element_id);
    Medium const* getMedium(std::size_t element_id) const;
    void checkElementHasMedium(std::size_t const element_id) const;

private:
    std::map<int, std::shared_ptr<Medium>> const& media_;
    MeshLib::PropertyVector<int> const* const material_ids_;
};
}  // namespace MaterialPropertyLib
