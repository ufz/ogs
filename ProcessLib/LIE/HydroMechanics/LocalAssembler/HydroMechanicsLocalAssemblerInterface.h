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

#include <utility>
#include <vector>

#include "BaseLib/Error.h"

#include "MeshLib/Elements/Element.h"

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
class HydroMechanicsLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    HydroMechanicsLocalAssemblerInterface(MeshLib::Element const& element,
                                          bool const is_axially_symmetric,
                                          std::size_t n_local_size,
                                          std::vector<unsigned>
                                              dofIndex_to_localIndex)
        : element_(element),
          is_axially_symmetric_(is_axially_symmetric),
          dofIndex_to_localIndex_(std::move(dofIndex_to_localIndex))
    {
        local_u_.resize(n_local_size);
        local_udot_.resize(n_local_size);
        local_b_.resize(local_u_.size());
        local_J_.resize(local_u_.size(), local_u_.size());
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "HydroMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x_,
                              std::vector<double> const& local_xdot_,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_dof_size = local_x_.size();

        local_u_.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            local_u_[dofIndex_to_localIndex_[i]] = local_x_[i];
        }
        local_udot_.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            local_udot_[dofIndex_to_localIndex_[i]] = local_xdot_[i];
        }
        local_b_.setZero();
        local_J_.setZero();

        assembleWithJacobianConcrete(t, dt, local_u_, local_udot_, local_b_,
                                     local_J_);

        local_b_data.resize(local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            local_b_data[i] = local_b_[dofIndex_to_localIndex_[i]];
        }

        local_Jac_data.resize(local_dof_size * local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            for (unsigned j = 0; j < local_dof_size; j++)
            {
                local_Jac_data[i * local_dof_size + j] = local_J_(
                    dofIndex_to_localIndex_[i], dofIndex_to_localIndex_[j]);
            }
        }
    }

    void postTimestepConcrete(std::vector<double> const& local_x_,
                              const double t, double const dt) override
    {
        auto const local_dof_size = local_x_.size();

        local_u_.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            local_u_[dofIndex_to_localIndex_[i]] = local_x_[i];
        }

        postTimestepConcreteWithVector(t, dt, local_u_);
    }

protected:
    virtual void assembleWithJacobianConcrete(double const t, double const dt,
                                              Eigen::VectorXd const& local_u,
                                              Eigen::VectorXd const& local_udot,
                                              Eigen::VectorXd& local_b,
                                              Eigen::MatrixXd& local_J) = 0;

    virtual void postTimestepConcreteWithVector(
        double const t, double const dt, Eigen::VectorXd const& local_u) = 0;

    MeshLib::Element const& element_;
    bool const is_axially_symmetric_;

private:
    Eigen::VectorXd local_u_;
    Eigen::VectorXd local_udot_;
    Eigen::VectorXd local_b_;
    Eigen::MatrixXd local_J_;
    // a vector for mapping the index in the local DoF vector to the index in
    // the complete local solution vector which also include nodes where DoF are
    // not assigned.
    std::vector<unsigned> const dofIndex_to_localIndex_;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
