/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "StokesFlowProcessData.h"

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"

namespace ProcessLib::StokesFlow
{
template <typename ShapeFunctionLiquidVelocity, typename ShapeFunctionPressure,
          typename IntegrationMethod, int GlobalDim>
class LocalAssemblerData : public StokesFlowLocalAssemblerInterface
{
    static const int liquid_velocity_index = 0;
    static const int pressure_index =
        ShapeFunctionLiquidVelocity::NPOINTS * GlobalDim;

    static const int liquid_velocity_size =
        ShapeFunctionLiquidVelocity::NPOINTS * GlobalDim;
    static const int pressure_size = ShapeFunctionPressure::NPOINTS;

    using ShapeMatrixTypeLiquidVelocity =
        ShapeMatrixPolicyType<ShapeFunctionLiquidVelocity, GlobalDim>;
    using ShapeMatrixTypePressure =
        ShapeMatrixPolicyType<ShapeFunctionPressure, GlobalDim>;

    using NodalVectorType = typename ShapeMatrixTypePressure::NodalVectorType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const /*local_matrix_size*/,
                       bool const is_axially_symmetric,
                       unsigned const integration_order,
                       StokesFlowProcessData const& process_data)
        : _element(element),
          _is_axially_symmetric(is_axially_symmetric),
          _integration_method(integration_order),
          _process_data(process_data)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices_v =
            NumLib::initShapeMatrices<ShapeFunctionLiquidVelocity,
                                      ShapeMatrixTypeLiquidVelocity, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        auto const shape_matrices_p =
            NumLib::initShapeMatrices<ShapeFunctionPressure,
                                      ShapeMatrixTypePressure, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices_v[ip].N, shape_matrices_v[ip].dNdx,
                shape_matrices_p[ip].N, shape_matrices_p[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices_v[ip].integralMeasure *
                    shape_matrices_v[ip].detJ);

            auto& ip_data = _ip_data[ip];

            ip_data.N_v_op = ShapeMatrixTypeLiquidVelocity::template MatrixType<
                GlobalDim, liquid_velocity_size>::Zero(GlobalDim,
                                                       liquid_velocity_size);

            ip_data.dNdx_v_op =
                ShapeMatrixTypeLiquidVelocity::template MatrixType<
                    GlobalDim * GlobalDim,
                    liquid_velocity_size>::Zero(GlobalDim * GlobalDim,
                                                liquid_velocity_size);

            for (int i = 0; i < GlobalDim; ++i)
            {
                ip_data.N_v_op
                    .template block<1, liquid_velocity_size / GlobalDim>(
                        i, i * liquid_velocity_size / GlobalDim)
                    .noalias() = shape_matrices_v[ip].N;

                ip_data.dNdx_v_op
                    .template block<GlobalDim,
                                    liquid_velocity_size / GlobalDim>(
                        i * GlobalDim, i * liquid_velocity_size / GlobalDim)
                    .noalias() = shape_matrices_v[ip].dNdx;
            }
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = liquid_velocity_size + pressure_size;
        assert(local_x.size() == local_matrix_size);

        auto local_p = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);

        auto local_K = MathLib::createZeroedMatrix<
            typename ShapeMatrixTypeLiquidVelocity::template MatrixType<
                local_matrix_size, local_matrix_size>>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<
            typename ShapeMatrixTypeLiquidVelocity::template VectorType<
                local_matrix_size>>(local_b_data, local_matrix_size);

        // Get block matrices
        // Kvv
        auto Kvv =
            local_K.template block<liquid_velocity_size, liquid_velocity_size>(
                liquid_velocity_index, liquid_velocity_index);
        // Kvp
        auto Kvp = local_K.template block<liquid_velocity_size, pressure_size>(
            liquid_velocity_index, pressure_index);
        // kpv
        auto Kpv = local_K.template block<pressure_size, liquid_velocity_size>(
            pressure_index, liquid_velocity_index);
        // bv
        auto bv = local_b.template segment<liquid_velocity_size>(
            liquid_velocity_index);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialPropertyLib::VariableArray vars;

        // create unit vector
        assert(GlobalDim == 2);
        auto const identity_mat =
            Eigen::Matrix<double, GlobalDim, GlobalDim>::Identity().eval();
        auto const unit_vector = Eigen::Map<const Eigen::RowVectorXd>(
                                     identity_mat.data(), identity_mat.size())
                                     .eval();

        // Get material properties
        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto& ip_data = _ip_data[ip];

            auto const& N_p = ip_data.N_p;
            auto const& N_v_op = ip_data.N_v_op;

            auto const& dNdx_v = ip_data.dNdx_v;
            auto const& dNdx_v_op = ip_data.dNdx_v_op;
            auto const& dNdx_p = ip_data.dNdx_p;

            auto const& w = ip_data.integration_weight;

            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_p, N_p, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // Use the viscosity model to compute the viscosity
            auto const mu = phase[MaterialPropertyLib::PropertyType::viscosity]
                                .template value<double>(vars, pos, t, dt);

            // Kvv
            if (_process_data.use_stokes_brinkman_form)
            {
                // permeability
                auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    medium[MaterialPropertyLib::PropertyType::permeability]
                        .value(vars, pos, t, dt));

                Kvv.noalias() +=
                    w * N_v_op.transpose() * mu * K.inverse() * N_v_op;
            }

            for (int i = 0; i < GlobalDim; ++i)
            {
                Kvv.template block<liquid_velocity_size / GlobalDim,
                                   liquid_velocity_size / GlobalDim>(
                       i * liquid_velocity_size / GlobalDim,
                       i * liquid_velocity_size / GlobalDim)
                    .noalias() += w * dNdx_v.transpose() * mu * dNdx_v;
            }

            // Kvp
            Kvp.noalias() += w * N_v_op.transpose() * dNdx_p;

            // Kpv
            Kpv.noalias() += w * N_p.transpose() * unit_vector * dNdx_v_op;

            // bv
            bv.noalias() += w * N_v_op.transpose() * b;
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_p = _ip_data[integration_point].N_p;

        // assumes N_p is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_p.data(), N_p.size());
    }

    void computeSecondaryVariableConcrete(
        double const /*t*/,
        double const /*dt*/,
        Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& /*local_x_dot*/) override
    {
        auto const local_p =
            local_x.template segment<pressure_size>(pressure_index);

        NumLib::interpolateToHigherOrderNodes<
            ShapeFunctionPressure,
            typename ShapeFunctionLiquidVelocity::MeshElement, GlobalDim>(
            _element, _is_axially_symmetric, local_p,
            *_process_data.pressure_interpolated);
    }

private:
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;
    IntegrationMethod const _integration_method;
    StokesFlowProcessData const& _process_data;

    std::vector<IntegrationPointData<ShapeMatrixTypeLiquidVelocity,
                                     ShapeMatrixTypePressure, GlobalDim,
                                     ShapeFunctionLiquidVelocity::NPOINTS>,
                Eigen::aligned_allocator<IntegrationPointData<
                    ShapeMatrixTypeLiquidVelocity, ShapeMatrixTypePressure,
                    GlobalDim, ShapeFunctionLiquidVelocity::NPOINTS>>>
        _ip_data;
};

}  // namespace ProcessLib::StokesFlow
