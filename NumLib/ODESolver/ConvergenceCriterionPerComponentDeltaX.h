/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/LinAlgEnums.h"
#include "ConvergenceCriterionPerComponent.h"

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
        const MathLib::VecNormType norm_type);

    bool hasDeltaXCheck() const override { return true; }
    bool hasResidualCheck() const override { return false; }

    void checkDeltaX(const GlobalVector& minus_delta_x,
                     GlobalVector const& x) override;
    void checkResidual(const GlobalVector& /*residual*/) override {}

    void reset() override { this->_satisfied = true; }

    void setDOFTable(const LocalToGlobalIndexMap& dof_table,
                     MeshLib::Mesh const& mesh) override;

private:
    const std::vector<double> _abstols;
    const std::vector<double> _reltols;
    LocalToGlobalIndexMap const* _dof_table = nullptr;
    MeshLib::Mesh const* _mesh = nullptr;
};

std::unique_ptr<ConvergenceCriterionPerComponentDeltaX>
createConvergenceCriterionPerComponentDeltaX(BaseLib::ConfigTree const& config);

}  // namespace NumLib
