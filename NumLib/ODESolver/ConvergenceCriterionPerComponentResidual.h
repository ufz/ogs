/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
        const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return true; }

    /// The function will only do diagnostic output and no actual check of the
    /// solution increment is made
    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& residual) override;

    void setDOFTable(const LocalToGlobalIndexMap& dof_table,
                     MeshLib::Mesh const& mesh) override;

private:
    const std::vector<double> abstols_;
    const std::vector<double> reltols_;
    LocalToGlobalIndexMap const* dof_table_ = nullptr;
    MeshLib::Mesh const* mesh_ = nullptr;
    std::vector<double> residual_norms_0_;
};

std::unique_ptr<ConvergenceCriterionPerComponentResidual>
createConvergenceCriterionPerComponentResidual(
    BaseLib::ConfigTree const& config);

}  // namespace NumLib
