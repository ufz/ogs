/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ForwardDifferencesJacobianAssembler.h"

#include "BaseLib/Error.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
ForwardDifferencesJacobianAssembler::ForwardDifferencesJacobianAssembler(
    std::vector<double>&& absolute_epsilons)
    : _absolute_epsilons(std::move(absolute_epsilons))
{
    if (_absolute_epsilons.empty())
    {
        OGS_FATAL("No values for the absolute epsilons have been given.");
    }
}

void ForwardDifferencesJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, const double t, double const dt,
    const std::vector<double>& local_x_data,
    const std::vector<double>& local_xdot_data,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    // TODO do not check in every call.
    if (local_x_data.size() % _absolute_epsilons.size() != 0)
    {
        OGS_FATAL(
            "The number of specified epsilons ({:d}) and the number of local "
            "d.o.f.s ({:d}) do not match, i.e., the latter is not divisible by "
            "the former.",
            _absolute_epsilons.size(), local_x_data.size());
    }

    auto const num_r_c =
        static_cast<Eigen::MatrixXd::Index>(local_x_data.size());

    auto const local_x =
        MathLib::toVector<Eigen::VectorXd>(local_x_data, num_r_c);
    auto const local_xdot =
        MathLib::toVector<Eigen::VectorXd>(local_xdot_data, num_r_c);

    auto local_Jac =
        MathLib::createZeroedMatrix(local_Jac_data, num_r_c, num_r_c);

    auto const num_dofs_per_component =
        local_x_data.size() / _absolute_epsilons.size();

    // Assemble with unperturbed local x to get M0, K0, and b0 used in the
    // finite differences below.
    local_assembler.assemble(t, dt, local_x_data, local_xdot_data, local_M_data,
                             local_K_data, local_b_data);

    // Residual  res := M xdot + K x - b
    // Computing Jac := dres/dx
    //                = d(M xdot)/dx + d(K x)/dx - db/dx
    for (Eigen::MatrixXd::Index i = 0; i < num_r_c; ++i)
    {
        // assume that local_x_data is ordered by component.
        auto const component = i / num_dofs_per_component;
        auto const eps = _absolute_epsilons[component];

        // Assemble with perturbed local x.
        _local_x_perturbed_data = local_x_data;
        _local_x_perturbed_data[i] = local_x_data[i] + eps;
        auto const x_p = MathLib::toVector<Eigen::VectorXd>(
            _local_x_perturbed_data, num_r_c);

        _local_xdot_perturbed_data = local_xdot_data;
        _local_xdot_perturbed_data[i] = local_xdot_data[i] + eps / dt;
        auto const xdot_p = MathLib::toVector<Eigen::VectorXd>(
            _local_xdot_perturbed_data, num_r_c);

        local_assembler.assemble(t, dt, _local_x_perturbed_data,
                                 _local_xdot_perturbed_data, _local_M_data,
                                 _local_K_data, _local_b_data);

        if (!local_M_data.empty() && !_local_M_data.empty())
        {
            auto const local_M_0 =
                MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
            auto const local_M_p =
                MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                (local_M_p * xdot_p - local_M_0 * local_xdot) / eps;
            _local_M_data.clear();
        }
        if (!local_K_data.empty() && !_local_K_data.empty())
        {
            auto const local_K_0 =
                MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
            auto const local_K_p =
                MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                (local_K_p * x_p - local_K_0 * local_x) / eps;
            _local_K_data.clear();
        }
        if (!local_b_data.empty() && !_local_b_data.empty())
        {
            auto const local_b_0 =
                MathLib::toVector<Eigen::VectorXd>(local_b_data, num_r_c);
            auto const local_b_p =
                MathLib::toVector<Eigen::VectorXd>(_local_b_data, num_r_c);
            local_Jac.col(i).noalias() -= (local_b_p - local_b_0) / eps;
            _local_b_data.clear();
        }
    }

    // Assemble with unperturbed local x to keep the internal assembler state
    // for next iteration or time step.
    local_assembler.assemble(t, dt, local_x_data, local_xdot_data, local_M_data,
                             local_K_data, local_b_data);
}

}  // namespace ProcessLib
