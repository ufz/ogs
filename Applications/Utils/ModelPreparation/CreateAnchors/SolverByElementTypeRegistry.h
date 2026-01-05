// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <functional>
#include <memory>
#include <typeindex>
#include <unordered_map>

#include "BaseLib/Logging.h"
#include "NumLib/Fem/FiniteElement/ElementTraitsLagrange.h"

namespace ApplicationUtils
{
/**
 * Instantiates a given solver interface for all mesh element types (of
 * dimension >= 1) and provides the correct solver implementation for a given
 * mesh element.
 */
template <typename SolverInterface,
          template <typename ShapeFunction>
          class SolverImplementationTplTpl>
class SolverByElementTypeRegistry final
{
    using SolverIntfPtr = std::unique_ptr<SolverInterface>;

    static std::unordered_map<std::type_index, SolverIntfPtr> initSolvers()
    {
        std::unordered_map<std::type_index, SolverIntfPtr> solvers;

        auto at_least_1D = []<typename ET>(ET*)
        { return ET::Element::dimension >= 1; };

        auto create_and_insert_solver = [&solvers]<typename ET>(ET*)
        {
            using SolverImplementation =
                SolverImplementationTplTpl<typename ET::ShapeFunction>;

            solvers[std::type_index(typeid(typename ET::Element))] =
                std::make_unique<SolverImplementation>();
        };

        using ElementTraitsFiltered =
            decltype(BaseLib::TMP::filter<NumLib::AllElementTraitsLagrange>(
                at_least_1D));

        BaseLib::TMP::foreach<ElementTraitsFiltered>(create_and_insert_solver);

        return solvers;
    }

    static const std::unordered_map<std::type_index, SolverIntfPtr> solvers_;

public:
    static SolverInterface const& getFor(MeshLib::Element const& e)
    {
        auto const type_idx = std::type_index(typeid(e));
        auto const it_type_solver = solvers_.find(type_idx);

        if (it_type_solver != solvers_.end())
        {
            return *it_type_solver->second;
        }
        OGS_FATAL(
            "You are trying to get a solver for an unknown mesh element "
            "type "
            "({:s}).",
            type_idx.name());
    }
};

template <typename SolverInterface,
          template <typename /* shp fct */>
          class SolverImplementationTplTpl>
const std::unordered_map<std::type_index, std::unique_ptr<SolverInterface>>
    SolverByElementTypeRegistry<SolverInterface,
                                SolverImplementationTplTpl>::solvers_{
        SolverByElementTypeRegistry<SolverInterface,
                                    SolverImplementationTplTpl>::initSolvers()};

}  // namespace ApplicationUtils
