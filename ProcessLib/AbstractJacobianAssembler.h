// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    virtual void checkPerturbationSize(
        int const max_non_deformation_dofs_per_node) const = 0;

    virtual void setNonDeformationComponentIDs(
        std::vector<int> const& non_deformation_component_ids) = 0;

    virtual void setNonDeformationComponentIDsNoSizeCheck(
        std::vector<int> const& non_deformation_component_ids) = 0;

    virtual bool needsPicardAssembly() const = 0;
};

}  // namespace ProcessLib
