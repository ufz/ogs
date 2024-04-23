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

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/Reflection/ReflectionIPData.h"

namespace ProcessLib
{
struct CellAverageData
{
    explicit CellAverageData(MeshLib::Mesh& mesh) : mesh_{mesh} {}

    template <typename LAIntf>
    void computeSecondaryVariable(
        int const dim,
        std::vector<std::unique_ptr<LAIntf>> const& local_assemblers)
    {
        auto const callback =
            [this, &local_assemblers](std::string const& name,
                                      unsigned const num_comp,
                                      auto&& flattened_ip_data_accessor)
        {
            computeCellAverages(name, num_comp, flattened_ip_data_accessor,
                                local_assemblers);
        };

        if (dim == 2)
        {
            ProcessLib::Reflection::
                forEachReflectedFlattenedIPDataAccessor<2, LAIntf>(
                    LAIntf::getReflectionDataForOutput(), callback);
        }
        else if (dim == 3)
        {
            ProcessLib::Reflection::
                forEachReflectedFlattenedIPDataAccessor<3, LAIntf>(
                    LAIntf::getReflectionDataForOutput(), callback);
        }
        else
        {
            OGS_FATAL(
                "Generic averaged output is only implemented for dimensions 2 "
                "and 3 ATM. You requested dim = {}.",
                dim);
        }
    }

private:
    MeshLib::PropertyVector<double>& getOrCreatePropertyVector(
        std::string const& name, unsigned const num_comp);

    void computeCellAverages(std::string const& name,
                             unsigned const num_comp,
                             auto&& flattened_ip_data_accessor,
                             auto const& local_assemblers)
    {
        auto& prop_vec = getOrCreatePropertyVector(name, num_comp);

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

    MeshLib::Mesh const& mesh_;
    std::map<std::string, MeshLib::PropertyVector<double>*> cell_averages_;
};
}  // namespace ProcessLib
