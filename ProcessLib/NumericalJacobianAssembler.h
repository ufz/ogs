// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "AbstractJacobianAssembler.h"

namespace ProcessLib
{
//! Base class for numerical Jacobian assemblers
class NumericalJacobianAssembler : public AbstractJacobianAssembler
{
public:
    //! Constructs a new instance for the numerical Jacobian.
    //! \param absolute_epsilons perturbations of the components of the local
    //! solution vector used    for evaluating the finite differences.
    explicit NumericalJacobianAssembler(std::vector<double>&& absolute_epsilons)
        : absolute_epsilons_(std::move(absolute_epsilons))
    {
        if (absolute_epsilons_.empty())
        {
            OGS_FATAL(
                "Numerical Jacobian assembler requires perturbation values "
                "(epsilons) for finite difference approximation, but none were "
                "provided.\n"
                "Please specify <epsilons> in the <jacobian_assembler> section "
                "of your project file.\n"
                "Example: <epsilons>1e-8 1e-8</epsilons> for a process with 2 "
                "non-deformation variables.");
        }
    }

    bool isPerturbationEnabled() const override
    {
        return true;
        // return !absolute_epsilons_.empty();
    }

    void checkPerturbationSize(
        int const max_non_deformation_dofs_per_node) const override
    {
        int const num_abs_eps = static_cast<int>(absolute_epsilons_.size());
        if (num_abs_eps != max_non_deformation_dofs_per_node)
        {
            OGS_FATAL(
                "Mismatch in numerical Jacobian perturbation configuration:\n"
                "  Provided epsilons:  {:d}\n"
                "  Required epsilons:  {:d} (one per non-deformation variable "
                "component)\n\n"
                "Each non-deformation variable needs exactly one epsilon value "
                "for numerical differentiation.\n"
                "Deformation variables use analytical Jacobian and do not "
                "require epsilons.\n"
                "Please adjust the <epsilons> array in your project file to "
                "match the number of required components.",
                num_abs_eps, max_non_deformation_dofs_per_node);
        }
    }

    void setNonDeformationComponentIDs(
        std::vector<int> const& non_deformation_component_ids) override
    {
        // So far, the number of specified perturbations is checked inside this
        // function, assuming the TH2M and THM processes do not use the
        // staggered scheme. If the staggered scheme is supported in these
        // processes, call `checkPerturbationSize(...)` and then
        // `setNonDeformationComponentIDsNoSizeCheck(...)` instead.
        checkPerturbationSize(non_deformation_component_ids.size());

        non_deformation_component_ids_ = non_deformation_component_ids;
    }

    void setNonDeformationComponentIDsNoSizeCheck(
        std::vector<int> const& non_deformation_component_ids) override
    {
        non_deformation_component_ids_ = non_deformation_component_ids;
    }

protected:
    auto getVariableComponentEpsilonsView() const
    {
        assert(!non_deformation_component_ids_.empty());
        namespace rv = ranges::views;
        return non_deformation_component_ids_ |
               rv::transform(
                   [this](int comp_id) -> double const&
                   {
                       return absolute_epsilons_[static_cast<std::size_t>(
                           comp_id)];
                   });
    }

    /// Perturbations of the variable components used for evaluating finite
    /// differences (excluding deformation).
    /// \note The perturbations must be specified component-wise for each
    /// non-deformation variable. If the number of perturbations is not equal
    /// to the number of non-deformation variable components, a fatal error is
    /// issued. The deformation block of the local Jacobian matrix is calculated
    /// analytically; therefore, perturbations for deformation components are
    /// not required.
    std::vector<double> const absolute_epsilons_;

    /// IDs of the components that are not deformation variables. Used by the
    /// numerical Jacobian assembler. It is manipulated by processes that use
    /// the numerical Jacobian assembler. Therefore, it is thread-safe.
    std::vector<int> non_deformation_component_ids_;
};
}  // namespace ProcessLib
