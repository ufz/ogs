/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib::Deformation
{
template <typename LocalAssemblerInterface,
          typename AddSecondaryVariableCallback, int DisplacementDim>
void solidMaterialInternalToSecondaryVariables(
    std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<
                      DisplacementDim>>> const& solid_materials,
    AddSecondaryVariableCallback const& add_secondary_variable)
{
    // Collect the internal variables for all solid materials.
    std::vector<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::InternalVariable>
        internal_variables;
    for (auto const& material_id__solid_material : solid_materials)
    {
        auto const variables =
            material_id__solid_material.second->getInternalVariables();
        copy(begin(variables), end(variables),
             back_inserter(internal_variables));
    }

    // Register the internal variables.
    for (auto const& internal_variable : internal_variables)
    {
        auto const& name = internal_variable.name;
        auto const& fct = internal_variable.getter;
        auto const num_components = internal_variable.num_components;
        DBUG("Registering internal variable {:s}.", name);

        auto getIntPtValues =
            [fct, num_components](
                LocalAssemblerInterface const& loc_asm,
                const double /*t*/,
                std::vector<GlobalVector*> const& /*x*/,
                std::vector<
                    NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
                std::vector<double>& cache) -> std::vector<double> const& {
            const unsigned num_int_pts = loc_asm.getNumberOfIntegrationPoints();

            cache.clear();
            auto cache_mat = MathLib::createZeroedMatrix<Eigen::Matrix<
                double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
                cache, num_components, num_int_pts);

            // TODO avoid the heap allocation (one per finite element)
            std::vector<double> cache_column(num_int_pts);

            for (unsigned i = 0; i < num_int_pts; ++i)
            {
                auto const& state = loc_asm.getMaterialStateVariablesAt(i);

                auto const& int_pt_values = fct(state, cache_column);
                assert(int_pt_values.size() == num_components);
                auto const int_pt_values_vec = MathLib::toVector(int_pt_values);

                cache_mat.col(i).noalias() = int_pt_values_vec;
            }

            return cache;
        };

        add_secondary_variable(name, num_components, std::move(getIntPtValues));
    }

}

}  // namespace ProcessLib::Deformation
