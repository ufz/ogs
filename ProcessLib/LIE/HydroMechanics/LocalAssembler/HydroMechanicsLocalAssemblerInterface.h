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
#include "MeshLib/Elements/Element.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{

/// Used for the extrapolation of the integration point values.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

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
        : _element(element),
          _is_axially_symmetric(is_axially_symmetric),
          _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
    {
        _local_u.resize(n_local_size);
        _local_u_prev.resize(n_local_size);
        _local_b.resize(_local_u.size());
        _local_J.resize(_local_u.size(), _local_u.size());
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev*/,
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
                              std::vector<double> const& local_x_prev_,
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
        _local_u_prev.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
        {
            _local_u_prev[_dofIndex_to_localIndex[i]] = local_x_prev_[i];
        }
        _local_b.setZero();
        _local_J.setZero();

        assembleWithJacobianConcrete(t, dt, _local_u, _local_u_prev, _local_b,
                                     _local_J);

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

    void postTimestepConcrete(Eigen::VectorXd const& local_x_,
                              Eigen::VectorXd const& /*local_x_prev*/,
                              const double t, double const dt,
                              int const /*process_id*/) override
    {
        _local_u.setZero();
        for (Eigen::Index i = 0; i < local_x_.rows(); i++)
        {
            _local_u[_dofIndex_to_localIndex[i]] = local_x_[i];
        }

        postTimestepConcreteWithVector(t, dt, _local_u);
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

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtFractureVelocity(
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

    virtual std::vector<double> const& getIntPtFracturePermeability(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

protected:
    virtual void assembleWithJacobianConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_u,
        Eigen::VectorXd const& local_u_prev, Eigen::VectorXd& local_b,
        Eigen::MatrixXd& local_J) = 0;

    virtual void postTimestepConcreteWithVector(
        double const t, double const dt, Eigen::VectorXd const& local_u) = 0;

    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;

private:
    Eigen::VectorXd _local_u;
    Eigen::VectorXd _local_u_prev;
    Eigen::VectorXd _local_b;
    Eigen::MatrixXd _local_J;
    // a vector for mapping the index in the local DoF vector to the index in
    // the complete local solution vector which also include nodes where DoF are
    // not assigned.
    std::vector<unsigned> const _dofIndex_to_localIndex;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
