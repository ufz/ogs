/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
             std::vector<std::string>{},
             BaseLib::easyBind(&CachedSecondaryVariable::getValue, this)}};
}

double CachedSecondaryVariable::getValue() const
{
    if (_needs_recomputation)
        evalFieldNoArgs();
    return _cached_nodal_values.get(_context.index);
}

SecondaryVariableFunctions CachedSecondaryVariable::getExtrapolator()
{
    // TODO copied from makeExtrapolator()
    auto const eval_residuals = [this](
        const double t,
        GlobalVector const& x,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::unique_ptr<GlobalVector> & /*result_cache*/
        ) -> GlobalVector const& {
        _extrapolator.calculateResiduals(1, *_extrapolatables, t, x, dof_table);
        return _extrapolator.getElementResiduals();
    };
    return {1, BaseLib::easyBind(&CachedSecondaryVariable::evalField, this),
            eval_residuals};
}

GlobalVector const& CachedSecondaryVariable::evalField(
    const double t,
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::unique_ptr<GlobalVector>& /*result_cache*/
    ) const
{
    // _current_solution = &x;
    // _dof_table = &dof_table;
    return evalFieldNoArgs();
}

GlobalVector const& CachedSecondaryVariable::evalFieldNoArgs() const
{
    if (!_needs_recomputation) {
        DBUG("%s does not need to be recomputed. Returning cached values",
             _internal_variable_name.c_str());
        return _cached_nodal_values;
    }
    DBUG("Recomputing %s.", _internal_variable_name.c_str());
    // TODO fix
    const double t = 0.0;
    _extrapolator.extrapolate(
        1, *_extrapolatables, t, *_current_solution, *_dof_table);
    auto const& nodal_values = _extrapolator.getNodalValues();
    MathLib::LinAlg::copy(nodal_values, _cached_nodal_values);
    _needs_recomputation = false;
    return nodal_values;
}

}  // namespace ProcessLib
