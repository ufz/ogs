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

#include "ProcessLib/Output/SecondaryVariable.h"
#include "ReflectionIPData.h"

namespace ProcessLib::Reflection
{
/// Adds secondary variables for all IP data obtained recursively from the given
/// \c reflection_data to the given secondary variable collection.
template <int Dim, typename LocAsmIF, typename ReflData>
void addReflectedSecondaryVariables(
    ReflData const& reflection_data,
    SecondaryVariableCollection& secondary_variables,
    NumLib::Extrapolator& extrapolator,
    std::vector<std::unique_ptr<LocAsmIF>> const& local_assemblers)
{
    forEachReflectedFlattenedIPDataAccessor<Dim, LocAsmIF>(
        reflection_data,
        [&secondary_variables, &local_assemblers, &extrapolator](
            std::string const& name,
            unsigned const num_comp,
            auto&& flattened_ip_data_accessor)
        {
            secondary_variables.addSecondaryVariable(
                name,
                makeExtrapolator2(num_comp, extrapolator, local_assemblers,
                                  flattened_ip_data_accessor));
        });
}
}  // namespace ProcessLib::Reflection
