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

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
class SmallDeformationLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    SmallDeformationLocalAssemblerInterface() : dofIndex_to_localIndex_{} {}
    SmallDeformationLocalAssemblerInterface(
        std::size_t n_local_size, std::vector<unsigned> dofIndex_to_localIndex)
        : dofIndex_to_localIndex_(std::move(dofIndex_to_localIndex))
    {
        local_u_.resize(n_local_size);
        local_b_.resize(local_u_.size());
        local_J_.resize(local_u_.size(), local_u_.size());
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x_,
                              std::vector<double> const& /*local_xdot*/,
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
        local_b_.setZero();
        local_J_.setZero();

        assembleWithJacobian(t, dt, local_u_, local_b_, local_J_);

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

    virtual void assembleWithJacobian(double const /*t*/, double const /*dt*/,
                                      Eigen::VectorXd const& /*local_u*/,
                                      Eigen::VectorXd& /*local_b*/,
                                      Eigen::MatrixXd& /*local_J*/)
    {
        OGS_FATAL(
            "SmallDeformationLocalAssemblerInterface::assembleWithJacobian() "
            "is not implemented");
    }

    void computeSecondaryVariableConcrete(
        double const t, double const /*dt*/, std::vector<double> const& local_x,
        std::vector<double> const& /*local_x_dot*/) override
    {
        if (!dofIndex_to_localIndex_.empty())
        {
            local_u_.setZero();
            for (std::size_t i = 0; i < local_x.size(); i++)
            {
                local_u_[dofIndex_to_localIndex_[i]] = local_x[i];
            }
        }

        computeSecondaryVariableConcreteWithVector(t, local_u_);
    }

    virtual std::vector<double> const& getIntPtSigmaXX(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYY(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaZZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXY(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaXZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtSigmaYZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXX(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYY(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonZZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXY(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonXZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilonYZ(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

protected:
    virtual void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_u) = 0;

private:
    Eigen::VectorXd local_u_;
    Eigen::VectorXd local_b_;
    Eigen::MatrixXd local_J_;
    std::vector<unsigned> const dofIndex_to_localIndex_;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
