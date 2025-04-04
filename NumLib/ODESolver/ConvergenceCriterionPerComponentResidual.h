/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "ConvergenceCriterionPerComponent.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

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
        std::vector<double>&& damping_alpha,
        bool daming_alpha_switch,
        const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return true; }
    bool hasNonNegativeDamping() const override
    {
        return _damping_alpha_switch;
    }

    /// The function will only do diagnostic output and no actual check of the
    /// solution increment is made
    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& residual) override;
    double getDampingFactor(GlobalVector const& minus_delta_x,
                            GlobalVector const& x,
                            double damping_orig) override;

    void setDOFTable(const LocalToGlobalIndexMap& dof_table,
                     MeshLib::Mesh const& mesh) override;

private:
    const std::vector<double> _abstols;
    const std::vector<double> _reltols;
    LocalToGlobalIndexMap const* _dof_table = nullptr;
    MeshLib::Mesh const* _mesh = nullptr;
    std::vector<double> _residual_norms_0;
    const std::vector<double> _damping_alpha;
    bool _damping_alpha_switch;
};

std::unique_ptr<ConvergenceCriterionPerComponentResidual>
createConvergenceCriterionPerComponentResidual(
    BaseLib::ConfigTree const& config);

}  // namespace NumLib
