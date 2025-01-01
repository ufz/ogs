/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "CellAverageData.h"
#include "ProcessLib/Reflection/ReflectionIPData.h"

namespace ProcessLib
{
namespace detail
{
void computeCellAverages(CellAverageData& cell_average_data,
                         std::string const& name, unsigned const num_comp,
                         auto&& flattened_ip_data_accessor,
                         auto const& local_assemblers)
{
    auto& prop_vec =
        cell_average_data.getOrCreatePropertyVector(name, num_comp);

    for (std::size_t i = 0; i < local_assemblers.size(); ++i)
    {
        auto const& loc_asm = *local_assemblers[i];
        auto const& ip_data = flattened_ip_data_accessor(loc_asm);
        assert(ip_data.size() % num_comp == 0);
        auto const num_ips =
            static_cast<Eigen::Index>(ip_data.size() / num_comp);
        Eigen::Map<const Eigen::MatrixXd> ip_data_mapped{ip_data.data(),
                                                         num_comp, num_ips};

        Eigen::Map<Eigen::VectorXd>{&prop_vec[i * num_comp], num_comp} =
            ip_data_mapped.rowwise().mean();
    }
}
}  // namespace detail

template <int dim, typename LAIntf>
void computeCellAverages(
    CellAverageData& cell_average_data,
    std::vector<std::unique_ptr<LAIntf>> const& local_assemblers)
{
    auto const callback = [&cell_average_data, &local_assemblers](
                              std::string const& name,
                              unsigned const num_comp,
                              auto&& flattened_ip_data_accessor)
    {
        detail::computeCellAverages(cell_average_data, name, num_comp,
                                    flattened_ip_data_accessor,
                                    local_assemblers);
    };

    ProcessLib::Reflection::forEachReflectedFlattenedIPDataAccessor<dim,
                                                                    LAIntf>(
        LAIntf::getReflectionDataForOutput(), callback);
}
}  // namespace ProcessLib
