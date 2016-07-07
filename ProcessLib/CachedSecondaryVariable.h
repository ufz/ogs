/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CACHEDSECONDARYVARIABLE_H
#define PROCESSLIB_CACHEDSECONDARYVARIABLE_H

#include "BaseLib/Functional.h"
#include "NumLib/Extrapolation/ExtrapolatableElementCollection.h"
#include "NumLib/NamedFunctionProvider.h"
#include "NumLib/NumericsConfig.h"
#include "SecondaryVariable.h"
#include "SecondaryVariableContext.h"

namespace ProcessLib
{
class CachedSecondaryVariable : public NumLib::NamedFunctionProvider
{
public:
    template <typename LocalAssemblerCollection,
              typename IntegrationPointValuesMethod>
    CachedSecondaryVariable(
        std::string const& internal_variable_name,
        NumLib::Extrapolator& extrapolator,
        LocalAssemblerCollection const& local_assemblers,
        IntegrationPointValuesMethod integration_point_values_method,
        SecondaryVariableContext const& context)
        : _extrapolator(extrapolator)
        , _extrapolatables(new NumLib::ExtrapolatableLocalAssemblerCollection<
                           LocalAssemblerCollection>{
              local_assemblers, integration_point_values_method})
        , _context(context)
        , _internal_variable_name(internal_variable_name)
    {
    }

    CachedSecondaryVariable(CachedSecondaryVariable const&) = delete;
    CachedSecondaryVariable(CachedSecondaryVariable&&) = delete;

    std::vector<NumLib::NamedFunction> getNamedFunctions() const override;

    SecondaryVariableFunctions getExtrapolator();

    void expire() { _needs_recomputation = true; }
private:
    double getValue() const;

    GlobalVector const& evalField(
        GlobalVector const& /*x*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::unique_ptr<GlobalVector>& /*result_cache*/
        ) const;

    GlobalVector const& evalFieldNoArgs() const;

    mutable GlobalVector _solid_density;
    mutable bool _needs_recomputation = true;

    NumLib::Extrapolator& _extrapolator;
    std::unique_ptr<NumLib::ExtrapolatableElementCollection> _extrapolatables;
    SecondaryVariableContext const& _context;
    std::string const _internal_variable_name;
};

}  // namespace ProcessLib

#endif // PROCESSLIB_CACHEDSECONDARYVARIABLE_H
