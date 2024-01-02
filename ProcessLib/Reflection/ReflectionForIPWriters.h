/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "ReflectionIPData.h"

namespace ProcessLib::Reflection
{
/// Adds IP data writers for all IP data obtained recursively from the given
/// \c reflection_data to the given IP writer vector.
template <int Dim, typename LocAsmIF, typename ReflData>
void addReflectedIntegrationPointWriters(
    ReflData const& reflection_data,
    std::vector<std::unique_ptr<MeshLib::IntegrationPointWriter>>&
        integration_point_writers,
    unsigned const integration_order,
    std::vector<std::unique_ptr<LocAsmIF>> const& local_assemblers)
{
    forEachReflectedFlattenedIPDataAccessor<Dim, LocAsmIF>(
        reflection_data,
        [&integration_point_writers, integration_order, &local_assemblers](
            std::string const& name,
            unsigned const num_comp,
            auto&& flattened_ip_data_accessor)
        {
            // TODO check if writer with such a name already exists.
            integration_point_writers.emplace_back(
                std::make_unique<MeshLib::IntegrationPointWriter>(
                    name + "_ip", num_comp, integration_order, local_assemblers,
                    flattened_ip_data_accessor));
        });
}
}  // namespace ProcessLib::Reflection
