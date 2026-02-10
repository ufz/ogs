// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "AbstractJacobianAssembler.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace ProcessLib
{
namespace detail
{
struct CompareJacobiansJacobianAssemblerImpl;
}  // namespace detail
//! Assembles the Jacobian matrix using two different Jacobian assemblers
//! and compares the assembled local Jacobian matrices.
//!
//! If the provided tolerances are exceeded, debugging information is logged in
//! the form of a Python script.
class CompareJacobiansJacobianAssembler final : public AbstractJacobianAssembler
{
    struct Key
    {
        // passkey idiom for "private" ctor and std::make_unique
    };

public:
    CompareJacobiansJacobianAssembler(
        std::unique_ptr<AbstractJacobianAssembler>&& asm1,
        std::unique_ptr<AbstractJacobianAssembler>&& asm2,
        double abs_tol,
        double rel_tol,
        bool fail_on_error,
        std::string const& log_file_path);

    explicit CompareJacobiansJacobianAssembler(
        std::shared_ptr<detail::CompareJacobiansJacobianAssemblerImpl> impl,
        Key);

    void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                              double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override;

    std::unique_ptr<AbstractJacobianAssembler> copy() const override;

    void checkPerturbationSize(
        int const max_non_deformation_dofs_per_node) const override;

    void setNonDeformationComponentIDs(
        std::vector<int> const& non_deformation_component_ids) override;

    void setNonDeformationComponentIDsNoSizeCheck(
        std::vector<int> const& non_deformation_component_ids) override;

    bool isPerturbationEnabled() const override;

private:
    // PIMPL idiom to enable copy()
    std::shared_ptr<detail::CompareJacobiansJacobianAssemblerImpl> impl_;
};

std::unique_ptr<CompareJacobiansJacobianAssembler>
createCompareJacobiansJacobianAssembler(BaseLib::ConfigTree const& config);

}  // namespace ProcessLib
