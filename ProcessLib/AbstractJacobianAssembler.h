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

#include <Eigen/Core>
#include <range/v3/view.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

#include "BaseLib/Error.h"

namespace ProcessLib
{
class LocalAssemblerInterface;

//! Base class for Jacobian assemblers.
class AbstractJacobianAssembler
{
public:
    //! Constructs a new instance for the numerical Jacobian.
    //! \param absolute_epsilons perturbations of the components of the local
    //! solution vector used    for evaluating the finite differences.
    explicit AbstractJacobianAssembler(
        std::vector<double> const&& absolute_epsilons)
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

    explicit AbstractJacobianAssembler() : absolute_epsilons_({}) {}

    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$.
    virtual void assembleWithJacobian(LocalAssemblerInterface& local_assembler,
                                      double const t, double const dt,
                                      std::vector<double> const& local_x,
                                      std::vector<double> const& local_x_prev,
                                      std::vector<double>& local_b_data,
                                      std::vector<double>& local_Jac_data) = 0;

    //! Assembles the Jacobian, the matrices \f$M\f$ and \f$K\f$, and the vector
    //! \f$b\f$ with coupling.
    virtual void assembleWithJacobianForStaggeredScheme(
        LocalAssemblerInterface& /*local_assembler*/, double const /*t*/,
        double const /*dt*/, Eigen::VectorXd const& /*local_x*/,
        Eigen::VectorXd const& /*local_x_prev*/, int const /*process_id*/,
        std::vector<double>& /*local_b_data*/,
        std::vector<double>& /*local_Jac_data*/)
    {
        // TODO make pure virtual.
        OGS_FATAL("not implemented.");
    }

    virtual std::unique_ptr<AbstractJacobianAssembler> copy() const = 0;

    virtual ~AbstractJacobianAssembler() = default;

    /**
     * Checks that the number of specified perturbations is not smaller
     * than the maximum number of non-deformation degrees of freedom per node.
     */
    void checkPerturbationSize(
        int const max_non_deformation_dofs_per_node) const
    {
        if (absolute_epsilons_.empty())
        {  // No epsilons specified — this assembler is used for the Picard
           // nonlinear solver or for the Newton nonlinear solver with
           // analytical Jacobian, so there is nothing to do.
            return;
        }

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
        std::vector<int> const& non_deformation_component_ids)
    {
        // Repeated in checkPerturbationSize to cover staggered scheme in TH2M
        // and THM, which may be supported in future.
        if (absolute_epsilons_.empty())
        {  // No epsilons specified — this assembler is used for the Picard
           // nonlinear solver or for the Newton nonlinear solver with
           // analytical Jacobian, so there is nothing to do.
            return;
        }

        // So far, the number of specified perturbations is checked inside this
        // function, assuming the TH2M and THM processes do not use the
        // staggered scheme. If the staggered scheme is supported in these
        // processes, call `checkPerturbationSize(...)` and then
        // `setNonDeformationComponentIDsNoSizeCheck(...)` instead.
        checkPerturbationSize(non_deformation_component_ids.size());

        non_deformation_component_ids_ = non_deformation_component_ids;
    }

    void setNonDeformationComponentIDsNoSizeCheck(
        std::vector<int> const& non_deformation_component_ids)
    {
        non_deformation_component_ids_ = non_deformation_component_ids;
    }

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

    auto isPerturbationEnabled() const { return !absolute_epsilons_.empty(); }

protected:
    /// IDs of the components that are not deformation variables. Used by the
    /// numerical Jacobian assembler. It is manipulated by processes that use
    /// the numerical Jacobian assembler. Therefore, it is thread-safe.
    std::vector<int> non_deformation_component_ids_;

    /// Perturbations of the variable components used for evaluating finite
    /// differences (excluding deformation).
    /// \note The perturbations must be specified component-wise for each
    /// non-deformation variable. If the number of perturbations is not equal
    /// to the number of non-deformation variable components, a fatal error is
    /// issued. The deformation block of the local Jacobian matrix is calculated
    /// analytically; therefore, perturbations for deformation components are
    /// not required.
    std::vector<double> const absolute_epsilons_;
};

}  // namespace ProcessLib
