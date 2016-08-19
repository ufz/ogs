/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_CONVERGENCECRITERIONPERCOMPONENTDELTAX_H
#define NUMLIB_CONVERGENCECRITERIONPERCOMPONENTDELTAX_H

#include "MathLib/LinAlg/LinAlgEnums.h"
#include "ConvergenceCriterionPerComponent.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

class ConvergenceCriterionPerComponentDeltaX
    : public ConvergenceCriterionPerComponent
{
public:
    explicit ConvergenceCriterionPerComponentDeltaX(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }

    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& /*residual*/) override {}

    void reset() override { _satisfied = true; }
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
};

std::unique_ptr<ConvergenceCriterionPerComponentDeltaX>
createConvergenceCriterionPerComponentDeltaX(BaseLib::ConfigTree const& config);

}  // namespace NumLib

#endif  // NUMLIB_CONVERGENCECRITERIONPERCOMPONENTDELTAX_H
