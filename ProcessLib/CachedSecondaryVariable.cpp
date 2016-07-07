/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CachedSecondaryVariable.h"
#include "BaseLib/Functional.h"
#include "MathLib/LinAlg/LinAlg.h"

namespace ProcessLib
{
std::vector<NumLib::NamedFunction> CachedSecondaryVariable::getNamedFunctions()
    const
{
    return {{_internal_variable_name,
             {},
             BaseLib::easyBind(&CachedSecondaryVariable::getValue, this)}};
}

double CachedSecondaryVariable::getValue() const
{
    if (_needs_recomputation)
        evalFieldNoArgs();
    return _solid_density.get(_context.getIndex());
}

SecondaryVariableFunctions CachedSecondaryVariable::getExtrapolator()
{
    // TODO copied from makeExtrapolator()
    auto const eval_residuals = [this](
        GlobalVector const& /*x*/,
        NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
        std::unique_ptr<GlobalVector> & /*result_cache*/
        ) -> GlobalVector const& {
        _extrapolator.calculateResiduals(*_extrapolatables);
        return _extrapolator.getElementResiduals();
    };
    return {BaseLib::easyBind(&CachedSecondaryVariable::evalField, this),
            eval_residuals};
}

GlobalVector const& CachedSecondaryVariable::evalField(
    GlobalVector const& /*x*/,
    NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
    std::unique_ptr<GlobalVector>& /*result_cache*/
    ) const
{
    return evalFieldNoArgs();
}

GlobalVector const& CachedSecondaryVariable::evalFieldNoArgs() const
{
    if (!_needs_recomputation) {
        INFO(
            "Solid density does not need to be recomputed. Returning "
            "cached values");
        return _solid_density;
    }
    DBUG("Recomputing solid density.");
    _extrapolator.extrapolate(*_extrapolatables);
    auto const& nodal_values = _extrapolator.getNodalValues();
    MathLib::LinAlg::copy(nodal_values, _solid_density);
    _needs_recomputation = false;
    return nodal_values;
}

}  // namespace ProcessLib
