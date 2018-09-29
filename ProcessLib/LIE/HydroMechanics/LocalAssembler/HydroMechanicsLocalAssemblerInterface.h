/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
        : _element(element),
          _is_axially_symmetric(is_axially_symmetric),
          _dofIndex_to_localIndex(std::move(dofIndex_to_localIndex))
    {
        _local_u.resize(n_local_size);
        _local_udot.resize(n_local_size);
        _local_b.resize(_local_u.size());
        _local_J.resize(_local_u.size(), _local_u.size());
    }

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_rhs_data*/) override
    {
        OGS_FATAL(
            "HydroMechanicsLocalAssembler: assembly without jacobian is not "
            "implemented.");
    }

    void assembleWithJacobian(double const t,
                              std::vector<double> const& local_x_,
                              std::vector<double> const& local_xdot_,
                              const double /*dxdot_dx*/, const double /*dx_dx*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_dof_size = local_x_.size();

        _local_u.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
            _local_u[_dofIndex_to_localIndex[i]] = local_x_[i];
        _local_udot.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
            _local_udot[_dofIndex_to_localIndex[i]] = local_xdot_[i];
        _local_b.setZero();
        _local_J.setZero();

        assembleWithJacobianConcrete(t, _local_u, _local_udot, _local_b,
                                     _local_J);

        local_b_data.resize(local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
            local_b_data[i] = _local_b[_dofIndex_to_localIndex[i]];

        local_Jac_data.resize(local_dof_size * local_dof_size);
        for (unsigned i = 0; i < local_dof_size; i++)
            for (unsigned j = 0; j < local_dof_size; j++)
                local_Jac_data[i * local_dof_size + j] = _local_J(
                    _dofIndex_to_localIndex[i], _dofIndex_to_localIndex[j]);
    }

    void computeSecondaryVariableConcrete(
        const double t, std::vector<double> const& local_x_) override
    {
        auto const local_dof_size = local_x_.size();

        _local_u.setZero();
        for (unsigned i = 0; i < local_dof_size; i++)
            _local_u[_dofIndex_to_localIndex[i]] = local_x_[i];

        computeSecondaryVariableConcreteWithVector(t, _local_u);
    }

protected:
    virtual void assembleWithJacobianConcrete(double const t,
                                              Eigen::VectorXd const& local_u,
                                              Eigen::VectorXd const& local_udot,
                                              Eigen::VectorXd& local_b,
                                              Eigen::MatrixXd& local_J) = 0;

    virtual void computeSecondaryVariableConcreteWithVector(
        double const t, Eigen::VectorXd const& local_u) = 0;

    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;

private:
    Eigen::VectorXd _local_u;
    Eigen::VectorXd _local_udot;
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
