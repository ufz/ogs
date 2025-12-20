// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ConvergenceCriterionPerComponent.h"
#include "MathLib/LinAlg/LinAlgEnums.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

//! Convergence criterion applying absolute or relative tolerances individually
//! to each component of the whole solution increment vector.
//!
//! A residual check is not done.
//! If both an absolute and a relative tolerances are specified, at least one of
//! them has to be satisfied.
class ConvergenceCriterionPerComponentDeltaX
    : public ConvergenceCriterionPerComponent
{
public:
    ConvergenceCriterionPerComponentDeltaX(
        std::vector<double>&& absolute_tolerances,
        std::vector<double>&& relative_tolerances,
        std::vector<double>&& damping_alpha,
        bool daming_alpha_switch,
        const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }
    bool hasNonNegativeDamping() const override
    {
        return _damping_alpha_switch;
    }

    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& /*residual*/) override {}
    double getDampingFactor(const GlobalVector& minus_delta_x,
                            GlobalVector const& x,
                            double damping_orig) override;

    void reset() override { this->_satisfied = true; }

    void setDOFTable(const LocalToGlobalIndexMap& dof_table,
                     MeshLib::Mesh const& mesh) override;

private:
    const std::vector<double> _abstols;
    const std::vector<double> _reltols;
    const std::vector<double> _damping_alpha;
    bool _damping_alpha_switch;
    LocalToGlobalIndexMap const* _dof_table = nullptr;
    MeshLib::Mesh const* _mesh = nullptr;
};

std::unique_ptr<ConvergenceCriterionPerComponentDeltaX>
createConvergenceCriterionPerComponentDeltaX(BaseLib::ConfigTree const& config);

}  // namespace NumLib
