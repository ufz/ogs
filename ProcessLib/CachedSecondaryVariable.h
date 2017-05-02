/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <utility>

#include "NumLib/Extrapolation/ExtrapolatableElementCollection.h"
#include "NumLib/NamedFunctionProvider.h"
#include "NumLib/NumericsConfig.h"
#include "SecondaryVariable.h"
#include "SecondaryVariableContext.h"

namespace ProcessLib
{
/*! Secondary variable which is extrapolated from integration points to mesh
 * nodes; the resulting extrapolated values are cached for subsequent use.
 */
class CachedSecondaryVariable final : public NumLib::NamedFunctionProvider
{
public:
    /*! Constructs a new instance.
     *
     * \param internal_variable_name the variable's name
     * \param extrapolator extrapolates integration point values to nodal
     * values.
     * \param local_assemblers provide the integration point values
     * \param integration_point_values_method extracts the integration point
     * values from the \c local_assemblers
     * \param context needed s.t. this class can act as a NamedFunction
     */
    template <typename LocalAssemblerCollection,
              typename IntegrationPointValuesMethod>
    CachedSecondaryVariable(
        std::string internal_variable_name,
        NumLib::Extrapolator& extrapolator,
        LocalAssemblerCollection const& local_assemblers,
        IntegrationPointValuesMethod integration_point_values_method,
        SecondaryVariableContext const& context)
        : _extrapolator(extrapolator),
          _extrapolatables(new NumLib::ExtrapolatableLocalAssemblerCollection<
                           LocalAssemblerCollection>{
              local_assemblers, integration_point_values_method}),
          _context(context),
          _internal_variable_name(std::move(internal_variable_name))
    {
    }

    CachedSecondaryVariable(CachedSecondaryVariable const&) = delete;
    CachedSecondaryVariable(CachedSecondaryVariable&&) = delete;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    //! Returns extrapolation functions that compute the secondary variable.
    SecondaryVariableFunctions getExtrapolator();

    //! Set that recomputation is necessary.
    void expire() { _needs_recomputation = true; }

private:
    //! Provides the value at the current index of the _context.
    double getValue() const;

    //! Computes the secondary Variable.
    GlobalVector const& evalField(
        GlobalVector const& /*x*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::unique_ptr<GlobalVector>& /*result_cache*/
        ) const;

    //! Computes the secondary Variable.
    GlobalVector const& evalFieldNoArgs() const;

    //! Cache for the computed values.
    mutable GlobalVector _cached_nodal_values;
    mutable bool _needs_recomputation = true;

    NumLib::Extrapolator& _extrapolator;
    std::unique_ptr<NumLib::ExtrapolatableElementCollection> _extrapolatables;
    SecondaryVariableContext const& _context;
    std::string const _internal_variable_name;
};

}  // namespace ProcessLib
