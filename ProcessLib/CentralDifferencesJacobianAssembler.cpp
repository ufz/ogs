/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CentralDifferencesJacobianAssembler.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace ProcessLib
{
CentralDifferencesJacobianAssembler::CentralDifferencesJacobianAssembler(
    std::vector<double>&& absolute_epsilons)
    : _absolute_epsilons(std::move(absolute_epsilons))
{
    if (_absolute_epsilons.empty())
    {
        OGS_FATAL("No values for the absolute epsilons have been given.");
    }
}

void CentralDifferencesJacobianAssembler::assembleWithJacobian(
    LocalAssemblerInterface& local_assembler, const double t, double const dt,
    const std::vector<double>& local_x_data,
    const std::vector<double>& local_x_prev_data,
    std::vector<double>& local_b_data, std::vector<double>& local_Jac_data)
{
    std::vector<double> local_M_data(local_Jac_data.size());
    std::vector<double> local_K_data(local_Jac_data.size());

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
    auto const local_x_prev =
        MathLib::toVector<Eigen::VectorXd>(local_x_prev_data, num_r_c);
    Eigen::VectorXd const local_xdot = (local_x - local_x_prev) / dt;

    auto local_Jac =
        MathLib::createZeroedMatrix(local_Jac_data, num_r_c, num_r_c);
    _local_x_perturbed_data = local_x_data;

    auto const num_dofs_per_component =
        local_x_data.size() / _absolute_epsilons.size();

    auto const vds = local_assembler.getVectorDeformationSegment();

    // Residual  res := M xdot + K x - b
    // Computing Jac := dres/dx
    //                = M dxdot/dx + dM/dx xdot + K dx/dx + dK/dx x - db/dx
    //                  with dxdot/dx = 1/dt and dx/dx = 1
    //                  (Note: dM/dx and dK/dx actually have the second and
    //                  third index transposed.)
    // The loop computes the dM/dx, dK/dx and db/dx terms, the rest is computed
    // afterwards.
    for (Eigen::MatrixXd::Index i = 0; i < num_r_c; ++i)
    {
        if ((vds != std::nullopt) && i >= vds->start_index &&
            (i < (vds->start_index + vds->size)))
        {
            // Avoid to compute the analytic block
            continue;
        }

        // assume that local_x_data is ordered by component.
        auto const component = i / num_dofs_per_component;
        auto const eps = _absolute_epsilons[component];

        _local_x_perturbed_data[i] += eps;
        local_assembler.assemble(t, dt, _local_x_perturbed_data,
                                 local_x_prev_data, local_M_data, local_K_data,
                                 local_b_data);

        _local_x_perturbed_data[i] = local_x_data[i] - eps;
        local_assembler.assemble(t, dt, _local_x_perturbed_data,
                                 local_x_prev_data, _local_M_data,
                                 _local_K_data, _local_b_data);

        _local_x_perturbed_data[i] = local_x_data[i];

        if (!local_M_data.empty())
        {
            auto const local_M_p =
                MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
            auto const local_M_m =
                MathLib::toMatrix(_local_M_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dM/dxi * x_dot
                (local_M_p - local_M_m) * local_xdot / (2.0 * eps);
            local_M_data.clear();
            _local_M_data.clear();
        }
        if (!local_K_data.empty())
        {
            auto const local_K_p =
                MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
            auto const local_K_m =
                MathLib::toMatrix(_local_K_data, num_r_c, num_r_c);
            local_Jac.col(i).noalias() +=
                // dK/dxi * x
                (local_K_p - local_K_m) * local_x / (2.0 * eps);
            local_K_data.clear();
            _local_K_data.clear();
        }
        if (!local_b_data.empty())
        {
            auto const local_b_p =
                MathLib::toVector<Eigen::VectorXd>(local_b_data, num_r_c);
            auto const local_b_m =
                MathLib::toVector<Eigen::VectorXd>(_local_b_data, num_r_c);
            local_Jac.col(i).noalias() -=
                // db/dxi
                (local_b_p - local_b_m) / (2.0 * eps);
            local_b_data.clear();
            _local_b_data.clear();
        }
    }

    // Assemble with unperturbed local x, i.e. compute M dxdot/dx + K dx/dx =
    // M/dt + K
    local_assembler.assemble(t, dt, local_x_data, local_x_prev_data,
                             local_M_data, local_K_data, local_b_data);

    // Compute remaining terms of the Jacobian.
    if (!local_M_data.empty())
    {
        auto local_M = MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_M / dt;
    }
    if (!local_K_data.empty())
    {
        auto local_K = MathLib::toMatrix(local_K_data, num_r_c, num_r_c);
        local_Jac.noalias() += local_K;
    }

    // Move the M and K contributions to the residuum for evaluation of nodal
    // forces, flow rates, and the like. Cleaning up the M's and K's storage so
    // it is not accounted for twice.
    auto b = [&]()
    {
        if (!local_b_data.empty())
        {
            return MathLib::toVector<Eigen::VectorXd>(local_b_data, num_r_c);
        }
        return MathLib::createZeroedVector<Eigen::VectorXd>(local_b_data,
                                                            num_r_c);
    }();

    if (!local_M_data.empty())
    {
        auto M = MathLib::toMatrix(local_M_data, num_r_c, num_r_c);
        b -= M * local_xdot;
        local_M_data.clear();
    }
    if (!local_K_data.empty())
    {
        auto K = MathLib::toMatrix(local_K_data, num_r_c, num_r_c);

        // Note: The deformation segment of \c b is already computed as
        // int{B^T sigma}dA, which is identical to K_uu * u. Therefore the
        // corresponding K block is set to zero.
        if (vds != std::nullopt)
        {
            K.block(vds->start_index, vds->start_index, vds->size, vds->size)
                .setZero();
        }

        b -= K * local_x;
        local_K_data.clear();
    }
}

std::unique_ptr<CentralDifferencesJacobianAssembler>
createCentralDifferencesJacobianAssembler(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__type}
    config.checkConfigParameter("type", "CentralDifferences");

    // TODO make non-optional.
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CentralDifferences__relative_epsilons}
    auto rel_eps = config.getConfigParameterOptional<std::vector<double>>(
        "relative_epsilons");
    //! \ogs_file_param{prj__processes__process__jacobian_assembler__CentralDifferences__component_magnitudes}
    auto comp_mag = config.getConfigParameterOptional<std::vector<double>>(
        "component_magnitudes");

    if (!!rel_eps != !!comp_mag)
    {
        OGS_FATAL(
            "Either both or none of <relative_epsilons> and "
            "<component_magnitudes> have to be specified.");
    }

    std::vector<double> abs_eps;

    if (rel_eps)
    {
        if (rel_eps->size() != comp_mag->size())
        {
            OGS_FATAL(
                "The numbers of components of  <relative_epsilons> and "
                "<component_magnitudes> do not match.");
        }

        abs_eps.resize(rel_eps->size());
        for (std::size_t i = 0; i < rel_eps->size(); ++i)
        {
            abs_eps[i] = (*rel_eps)[i] * (*comp_mag)[i];
        }
    }
    else
    {
        // By default 1e-8 is used as epsilon for all components.
        // TODO: remove this default value.
        abs_eps.emplace_back(1e-8);
    }

    return std::make_unique<CentralDifferencesJacobianAssembler>(
        std::move(abs_eps));
}

std::unique_ptr<AbstractJacobianAssembler>
CentralDifferencesJacobianAssembler::copy() const
{
    return std::make_unique<CentralDifferencesJacobianAssembler>(*this);
}

}  // namespace ProcessLib
