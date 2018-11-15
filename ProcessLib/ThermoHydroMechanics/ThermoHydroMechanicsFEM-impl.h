/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ThermoHydroMechanicsFEM.h"

#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include <iostream>

namespace ProcessLib
{
namespace ThermoHydroMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                   ShapeFunctionPressure, IntegrationMethod,
                                   DisplacementDim>::
    ThermoHydroMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ThermoHydroMechanicsProcessData<DisplacementDim>& process_data)
    : _process_data(process_data),
      _integration_method(integration_order),
      _element(e),
      _is_axially_symmetric(is_axially_symmetric)
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _ip_data.reserve(n_integration_points);
    _secondary_data.N_u.resize(n_integration_points);

    auto const shape_matrices_u =
        initShapeMatrices<ShapeFunctionDisplacement,
                          ShapeMatricesTypeDisplacement, IntegrationMethod,
                          DisplacementDim>(e, is_axially_symmetric,
                                           _integration_method);

    auto const shape_matrices_p =
        initShapeMatrices<ShapeFunctionPressure, ShapeMatricesTypePressure,
                          IntegrationMethod, DisplacementDim>(
            e, is_axially_symmetric, _integration_method);

    auto const& solid_material =
        MaterialLib::Solids::selectSolidConstitutiveRelation(
            _process_data.solid_materials, _process_data.material_ids,
            e.getID());

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        _ip_data.emplace_back(solid_material);
        auto& ip_data = _ip_data[ip];
        auto const& sm_u = shape_matrices_u[ip];
        _ip_data[ip].integration_weight =
            _integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        // Initialize current time step values
        static const int kelvin_vector_size =
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value;
        ip_data.sigma_eff.setZero(kelvin_vector_size);
        ip_data.eps.setZero(kelvin_vector_size);
        ip_data.eps_m.setZero(kelvin_vector_size);
        ip_data.eps_m_prev.resize(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        ip_data.eps_prev.resize(kelvin_vector_size);
        ip_data.sigma_eff_prev.resize(kelvin_vector_size);

        ip_data.N_u_op = ShapeMatricesTypeDisplacement::template MatrixType<
            DisplacementDim, displacement_size>::Zero(DisplacementDim,
                                                      displacement_size);
        for (int i = 0; i < DisplacementDim; ++i)
            ip_data.N_u_op
                .template block<1, displacement_size / DisplacementDim>(
                    i, i * displacement_size / DisplacementDim)
                .noalias() = sm_u.N;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        ip_data.N_p = shape_matrices_p[ip].N;
        ip_data.dNdx_p = shape_matrices_p[ip].dNdx;

        _secondary_data.N_u[ip] = shape_matrices_u[ip].N;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, DisplacementDim>::
    assembleWithJacobian(double const t, std::vector<double> const& local_x,
                         std::vector<double> const& local_xdot,
                         const double /*dxdot_dx*/, const double /*dx_dx*/,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    assert(local_x.size() ==
           pressure_size + displacement_size + temperature_size);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_index,
                                      displacement_size);

    auto T_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            temperature_size> const>(local_xdot.data() + temperature_index,
                                     temperature_size);

    auto p_dot =
        Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
            pressure_size> const>(local_xdot.data() + pressure_index,
                                  pressure_size);
    auto u_dot =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_xdot.data() + displacement_index,
                                      displacement_size);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            temperature_size + displacement_size + pressure_size,
            temperature_size + displacement_size + pressure_size>>(
        local_Jac_data, displacement_size + pressure_size + temperature_size,
        displacement_size + pressure_size + temperature_size);

    auto local_rhs = MathLib::createZeroedVector<
        typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size + pressure_size + temperature_size>>(
        local_rhs_data, displacement_size + pressure_size + temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType MTT;
    MTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTT_coeff;
    KTT_coeff.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTT;
    KTT.setZero(temperature_size, temperature_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTp;
    KTp.setZero(temperature_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType KTp_coeff;
    KTp_coeff.setZero(temperature_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType laplace_p;
    laplace_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_p;
    storage_p.setZero(pressure_size, pressure_size);

    typename ShapeMatricesTypePressure::NodalMatrixType storage_T;
    storage_T.setZero(pressure_size, temperature_size);

    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        Kup;
    Kup.setZero(displacement_size, pressure_size);
    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, temperature_size>
        KuT;
    KuT.setZero(displacement_size, temperature_size);

    double const& dt = _process_data.dt;

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& w = _ip_data[ip].integration_weight;

        auto const& N_u_op = _ip_data[ip].N_u_op;

        auto const& N_u = _ip_data[ip].N_u;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const& N_p = _ip_data[ip].N_p;
        auto const& dNdx_p = _ip_data[ip].dNdx_p;

        // same shape function for pressure and temperature since they have the
        // same order
        auto const& N_T = N_p;
        auto const& dNdx_T = dNdx_p;
        auto const T_int_pt = N_T * T;
        // auto const p_int_pt = N_T * p;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        auto& eps = _ip_data[ip].eps;
        auto const& sigma_eff = _ip_data[ip].sigma_eff;

        double const S = _process_data.specific_storage(t, x_position)[0];
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];
        double const alpha_s =
            _process_data.solid_linear_thermal_expansion_coefficient(
                t, x_position)[0];
        double const beta_f =
            _process_data.fluid_volumetric_thermal_expansion_coefficient(
                t, x_position)[0];
        double const lambda_f =
            _process_data.fluid_thermal_conductivity(t, x_position)[0];
        double const lambda_s =
            _process_data.solid_thermal_conductivity(t, x_position)[0];
        double const C_f =
            _process_data.fluid_specific_heat_capacity(t, x_position)[0];
        double const C_s =
            _process_data.solid_specific_heat_capacity(t, x_position)[0];
        double const T0 = _process_data.reference_temperature(t, x_position)[0];
        auto const alpha = _process_data.biot_coefficient(t, x_position)[0];
        auto const rho_sr = _process_data.solid_density(t, x_position)[0];
        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const porosity = _process_data.porosity(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;
        auto const& identity2 = MathLib::KelvinVector::Invariants<
            MathLib::KelvinVector::KelvinVectorDimensions<
                DisplacementDim>::value>::identity2;

        double const delta_T(T_int_pt - T0);
        double const thermal_strain = alpha_s * delta_T;

        double const rho_s = rho_sr * (1 - 3 * thermal_strain);

        auto velocity = (-K_over_mu * dNdx_p * p).eval();
        double const rho_f = rho_fr * (1 - beta_f * delta_T);
        velocity += K_over_mu * rho_f * b;

        //
        // displacement equation, displacement part
        //
        eps.noalias() = B * u;
        auto C = _ip_data[ip].updateConstitutiveRelationThermal(
            t, x_position, dt, u,
            _process_data.reference_temperature(t, x_position)[0],
            thermal_strain);

        local_Jac
            .template block<displacement_size, displacement_size>(
                displacement_index, displacement_index)
            .noalias() += B.transpose() * C * B * w;

        auto const rho = rho_s * (1 - porosity) + porosity * rho_f;
        local_rhs.template segment<displacement_size>(displacement_index)
            .noalias() -=
            (B.transpose() * sigma_eff - N_u_op.transpose() * rho * b) * w;

        //
        // displacement equation, pressure part (K_up)
        //
        Kup.noalias() += B.transpose() * alpha * identity2 * N_p * w;

        //
        // pressure equation, pressure part (K_pp and M_pp).
        //
        laplace_p.noalias() += dNdx_p.transpose() * K_over_mu * dNdx_p * w;

        storage_p.noalias() += N_p.transpose() * S * N_p * w;
        //
        //  RHS, pressure part
        //
        local_rhs.template segment<pressure_size>(pressure_index).noalias() +=
            dNdx_p.transpose() * rho_f * K_over_mu * b * w;
        //
        // pressure equation, temperature part (M_pT)
        //
        auto const beta = porosity * beta_f + (1 - porosity) * 3 * alpha_s;
        storage_T.noalias() += N_T.transpose() * beta * N_T * w;

        //
        // pressure equation, displacement part.
        //
        // Reusing Kup.transpose().

        //
        // temperature equation, temperature part.
        //
        auto const lambda = porosity * lambda_f + (1 - porosity) * lambda_s;
        KTT.noalias() += (dNdx_T.transpose() * lambda * dNdx_T +
                          dNdx_T.transpose() * velocity * N_p * rho_f * C_f) *
                         w;
        // coefficient matrix which is used for caculating the residual
        auto const heat_capacity =
            porosity * C_f * rho_f + (1 - porosity) * C_s * rho_sr;
        MTT.noalias() += N_T.transpose() * heat_capacity * N_T * w;

        //
        // temperature equation, pressure part
        //
        KTp.noalias() += K_over_mu * rho_f * C_f * N_T.transpose() *
                         (dNdx_T * T).transpose() * dNdx_T * w;
    }
    // temperature equation, temperature part
    local_Jac
        .template block<temperature_size, temperature_size>(temperature_index,
                                                            temperature_index)
        .noalias() += KTT + MTT / dt;

    // temperature equation, pressure part
    local_Jac
        .template block<temperature_size, pressure_size>(temperature_index,
                                                         pressure_index)
        .noalias() -= KTp;
    // displacement equation, temperature part
    local_Jac
        .template block<displacement_size, temperature_size>(displacement_index,
                                                             temperature_index)
        .noalias() -= KuT;

    // displacement equation, pressure part
    local_Jac
        .template block<displacement_size, pressure_size>(displacement_index,
                                                          pressure_index)
        .noalias() -= Kup;

    // pressure equation, temperature part.
    local_Jac
        .template block<pressure_size, temperature_size>(pressure_index,
                                                         temperature_index)
        .noalias() -= storage_T / dt;

    // pressure equation, pressure part.
    local_Jac
        .template block<pressure_size, pressure_size>(pressure_index,
                                                      pressure_index)
        .noalias() += laplace_p + storage_p / dt;

    // pressure equation, displacement part.
    local_Jac
        .template block<pressure_size, displacement_size>(pressure_index,
                                                          displacement_index)
        .noalias() += Kup.transpose() / dt;

    // pressure equation (f_p)
    local_rhs.template segment<pressure_size>(pressure_index).noalias() -=
        laplace_p * p + storage_p * p_dot - storage_T * T_dot +
        Kup.transpose() * u_dot;

    // displacement equation (f_u)
    local_rhs.template segment<displacement_size>(displacement_index)
        .noalias() += Kup * p;

    // temperature equation (f_T)
    local_rhs.template segment<temperature_size>(temperature_index).noalias() -=
        KTT * T + MTT * T_dot;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
std::vector<double> const& ThermoHydroMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunctionPressure, IntegrationMethod,
    DisplacementDim>::getIntPtDarcyVelocity(const double t,
                                            GlobalVector const&
                                                current_solution,
                                            NumLib::LocalToGlobalIndexMap const&
                                                dof_table,
                                            std::vector<double>& cache) const
{
    auto const num_intpts = _ip_data.size();

    auto const indices = NumLib::getIndices(_element.getID(), dof_table);
    assert(!indices.empty());
    auto const local_x = current_solution.get(indices);

    cache.clear();
    auto cache_matrix = MathLib::createZeroedMatrix<Eigen::Matrix<
        double, DisplacementDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, DisplacementDim, num_intpts);

    SpatialPosition pos;
    pos.setElementID(_element.getID());

    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    SpatialPosition x_position;
    x_position.setElementID(_element.getID());
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        double const K_over_mu =
            _process_data.intrinsic_permeability(t, x_position)[0] /
            _process_data.fluid_viscosity(t, x_position)[0];

        auto const rho_fr = _process_data.fluid_density(t, x_position)[0];
        auto const& b = _process_data.specific_body_force;

        // Compute the velocity
        auto const& dNdx_p = _ip_data[ip].dNdx_p;
        cache_matrix.col(ip).noalias() =
            -K_over_mu * dNdx_p * p - K_over_mu * rho_fr * b;
    }

    return cache;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, DisplacementDim>::
    postNonLinearSolverConcrete(std::vector<double> const& local_x,
                                double const t,
                                bool const use_monolithic_scheme)
{
    const int displacement_offset =
        use_monolithic_scheme ? displacement_index : 0;

    auto u =
        Eigen::Map<typename ShapeMatricesTypeDisplacement::template VectorType<
            displacement_size> const>(local_x.data() + displacement_offset,
                                      displacement_size);

    auto T = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        temperature_size> const>(local_x.data() + temperature_index,
                                 temperature_size);

    double const& dt = _process_data.dt;
    SpatialPosition x_position;
    x_position.setElementID(_element.getID());

    int const n_integration_points = _integration_method.getNumberOfPoints();
    for (int ip = 0; ip < n_integration_points; ip++)
    {
        x_position.setIntegrationPoint(ip);
        auto const& N_u = _ip_data[ip].N_u;
        auto const& N_T = _ip_data[ip].N_p;
        auto const& dNdx_u = _ip_data[ip].dNdx_u;

        auto const x_coord =
            interpolateXCoordinate<ShapeFunctionDisplacement,
                                   ShapeMatricesTypeDisplacement>(_element,
                                                                  N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, _is_axially_symmetric);

        double const T0 = _process_data.reference_temperature(t, x_position)[0];
        double const alpha_s =
            _process_data.solid_linear_thermal_expansion_coefficient(
                t, x_position)[0];

        double const T_int_pt = N_T * T;

        double const delta_T(T_int_pt - T0);
        double const thermal_strain = alpha_s * delta_T;

        auto& eps = _ip_data[ip].eps;
        eps.noalias() = B * u;

        _ip_data[ip].updateConstitutiveRelationThermal(
            t, x_position, dt, u,
            _process_data.reference_temperature(t, x_position)[0],
            thermal_strain);
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunctionPressure,
          typename IntegrationMethod, int DisplacementDim>
void ThermoHydroMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                        ShapeFunctionPressure,
                                        IntegrationMethod, DisplacementDim>::
    computeSecondaryVariableConcrete(double const /*t*/,
                                     std::vector<double> const& local_x)
{
    auto p = Eigen::Map<typename ShapeMatricesTypePressure::template VectorType<
        pressure_size> const>(local_x.data() + pressure_index, pressure_size);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunctionPressure, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(_element, _is_axially_symmetric, p,
                         *_process_data.pressure_interpolated);
}

}  // namespace ThermoHydroMechanics
}  // namespace ProcessLib
