/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>
#include "MathLib/LinAlg/LinAlgEnums.h"
#include "ConvergenceCriterionPerComponent.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

//! Convergence criterion applying absolute or relative tolerances individually
//! to each component of the whole residual vector.
//!
//! A check of the solution increment is not done.
//! If both an absolute and a relative tolerances are specified, at least one of
//! them has to be satisfied.
class ConvergenceCriterionPerComponentResidual
    : public ConvergenceCriterionPerComponent
{
public:
    ConvergenceCriterionPerComponentResidual(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return true; }

    /// The function will only do diagnostic output and no actual check of the
    /// solution increment is made
    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
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
