/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERIONPERCOMPONENTRESIDUAL_H
#define NUMLIB_CONVERGENCECRITERIONPERCOMPONENTRESIDUAL_H

#include <vector>
#include "MathLib/LinAlg/LinAlgEnums.h"
#include "ConvergenceCriterionPerComponent.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

class ConvergenceCriterionPerComponentResidual
    : public ConvergenceCriterionPerComponent
{
public:
    ConvergenceCriterionPerComponentResidual(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return false; }
    bool hasResidualCheck() const override { return true; }

    void checkDeltaX(const GlobalVector& /*minus_delta_x*/,
                     GlobalVector const& /*x*/) override {}
    void checkResidual(const GlobalVector& residual) override;

    void preFirstIteration() override { _is_first_iteration = true; }
    void reset() override { _satisfied = true; _is_first_iteration = false; }
    bool isSatisfied() const override { return _satisfied; }

    void setDOFTable(const LocalToGlobalIndexMap& dof_table,
                     MeshLib::Mesh const& mesh) override;

private:
    const std::vector<double> _abstols;
    const std::vector<double> _reltols;
    const MathLib::VecNormType _norm_type;
    bool _satisfied = true;
    LocalToGlobalIndexMap const* _dof_table = nullptr;
    MeshLib::Mesh const* _mesh = nullptr;
    bool _is_first_iteration = true;
    std::vector<double> _residual_norms_0;
};

std::unique_ptr<ConvergenceCriterionPerComponentResidual>
createConvergenceCriterionPerComponentResidual(
    BaseLib::ConfigTree const& config);

}  // namespace NumLib

#endif  // NUMLIB_CONVERGENCECRITERIONPERCOMPONENTRESIDUAL_H
