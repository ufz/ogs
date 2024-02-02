/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    SmallDeformationLocalAssemblerInterface() : _dofIndex_to_localIndex{} {}
    SmallDeformationLocalAssemblerInterface(
        std::size_t n_local_size, std::vector<unsigned> dofIndex_to_localIndex)
        : _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
    {
        _local_u.resize(n_local_size);
        _local_b.resize(_local_u.size());
        _local_J.resize(_local_u.size(), _local_u.size());
    }

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x_,
                              std::vector<double> const& /*local_x_prev*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_dof_size = local_x_.size();

        _local_u.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            _local_u[_dofIndex_to_localIndex[i]] = local_x_[i];
        }
        _local_b.setZero();
        _local_J.setZero();

        assembleWithJacobian(t, dt, _local_u, _local_b, _local_J);

        local_b_data.resize(local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            local_b_data[i] = _local_b[_dofIndex_to_localIndex[i]];
        }

        local_Jac_data.resize(local_dof_size * local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            for (unsigned j = 0; j < local_dof_size; j++)
            {
                local_Jac_data[i * local_dof_size + j] = _local_J(
                    _dofIndex_to_localIndex[i], _dofIndex_to_localIndex[j]);
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
        double const t, double const /*dt*/, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& /*local_x_prev*/) override
    {
        if (!_dofIndex_to_localIndex.empty())
        {
            _local_u.setZero();
            for (auto i = 0; i < local_x.rows(); i++)
            {
                _local_u[_dofIndex_to_localIndex[i]] = local_x[i];
            }
        }

        computeSecondaryVariableConcreteWithVector(t, _local_u);
    }

    virtual std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtFractureStress(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtFractureAperture(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

protected:
    virtual void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_u) = 0;

private:
    Eigen::VectorXd _local_u;
    Eigen::VectorXd _local_b;
    Eigen::MatrixXd _local_J;
    std::vector<unsigned> const _dofIndex_to_localIndex;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
