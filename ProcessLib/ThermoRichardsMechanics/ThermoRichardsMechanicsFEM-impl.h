/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on November 29, 2017, 2:03 PM
 */

#pragma once

#include <spdlog/fmt/bundled/format.h>

#include <cassert>
#include <type_traits>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MaterialLib/MPL/Utils/GetLiquidThermalExpansivity.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/Graph/Get.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidDensity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/LiquidViscosity.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "ProcessLib/Utils/TransposeInPlace.h"
#include "ThermoRichardsMechanicsFEM.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement, ShapeFunction,
                                      DisplacementDim, ConstitutiveTraits>::
    ThermoRichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
            process_data)
    : LocalAssemblerInterface<DisplacementDim, ConstitutiveTraits>(
          e, integration_method, is_axially_symmetric, process_data)
{
    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    ip_data_.resize(n_integration_points);

    auto const shape_matrices_u =
        NumLib::initShapeMatrices<ShapeFunctionDisplacement,
                                  ShapeMatricesTypeDisplacement,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                  DisplacementDim>(e, is_axially_symmetric,
                                                   integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& ip_data = ip_data_[ip];
        auto const& sm_u = shape_matrices_u[ip];
        ip_data_[ip].integration_weight =
            integration_method.getWeightedPoint(ip).getWeight() *
            sm_u.integralMeasure * sm_u.detJ;

        ip_data.N_u = sm_u.N;
        ip_data.dNdx_u = sm_u.dNdx;

        // ip_data.N_p and ip_data.dNdx_p are used for both p and T variables
        ip_data.N_p = shape_matrices[ip].N;
        ip_data.dNdx_p = shape_matrices[ip].dNdx;
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<
    ShapeFunctionDisplacement, ShapeFunction, DisplacementDim,
    ConstitutiveTraits>::setInitialConditionsConcrete(Eigen::VectorXd const
                                                          local_x,
                                                      double const t,
                                                      int const /*process_id*/)
{
    assert(local_x.size() ==
           temperature_size + pressure_size + displacement_size);

    auto const p_L = local_x.template segment<pressure_size>(pressure_index);
    auto const T =
        local_x.template segment<temperature_size>(temperature_index);

    constexpr double dt = std::numeric_limits<double>::quiet_NaN();
    auto const& medium =
        *this->process_data_.media_map.getMedium(this->element_.getID());

    MediaData const media_data{medium};

    typename ConstitutiveTraits::ConstitutiveSetting const constitutive_setting;
    typename ConstitutiveTraits::ConstitutiveModels models(
        this->process_data_, this->solid_material_);

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();
    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        // N is used for both T and p variables.
        auto const& N = ip_data_[ip].N_p;

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->element_.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data_[ip].N_u))};

        double T_ip;
        NumLib::shapeFunctionInterpolate(T, N, T_ip);

        double p_cap_ip;
        NumLib::shapeFunctionInterpolate(-p_L, N, p_cap_ip);

        MPL::VariableArray variables;
        variables.capillary_pressure = p_cap_ip;
        variables.liquid_phase_pressure = -p_cap_ip;
        // setting pG to 1 atm
        // TODO : rewrite equations s.t. p_L = pG-p_cap
        variables.gas_phase_pressure = 1.0e5;
        variables.temperature = T_ip;

        double const S_L =
            medium.property(MPL::PropertyType::saturation)
                .template value<double>(variables, x_position, t, dt);
        std::get<PrevState<SaturationData>>(this->prev_states_[ip])->S_L = S_L;

        constitutive_setting.init(models, t, dt, x_position, media_data,
                                  {T_ip, 0, {}}, this->current_states_[ip],
                                  this->prev_states_[ip]);

        if (this->process_data_.initial_stress.value)
        {
            variables.liquid_saturation = S_L;
            convertInitialStressType(ip, t, x_position, medium, variables,
                                     -p_cap_ip);
        }
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    convertInitialStressType(unsigned const ip,
                             double const t,
                             ParameterLib::SpatialPosition const x_position,
                             MaterialPropertyLib::Medium const& medium,
                             MPL::VariableArray const& variables,
                             double const p_at_ip)
{
    bool constexpr is_strain_temperature_constitutive =
        std::is_same<ConstitutiveStress_StrainTemperature::ConstitutiveTraits<
                         DisplacementDim>,
                     ConstitutiveTraits>::value;
    if (is_strain_temperature_constitutive &&
        this->process_data_.initial_stress.type ==
            InitialStress::Type::Effective)
    {
        return;
    }

    if (!is_strain_temperature_constitutive &&
        this->process_data_.initial_stress.type == InitialStress::Type::Total)
    {
        return;
    }

    double const alpha_b =
        medium.property(MPL::PropertyType::biot_coefficient)
            .template value<double>(variables, x_position, t, 0.0 /*dt*/);

    double const bishop =
        medium.property(MPL::PropertyType::bishops_effective_stress)
            .template value<double>(variables, x_position, t, 0.0 /*dt*/);

    ConstitutiveTraits::ConstitutiveSetting::convertInitialStressType(
        this->current_states_[ip], this->prev_states_[ip],
        bishop * alpha_b * p_at_ip * Invariants::identity2);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    assembleWithJacobian(double const t, double const dt,
                         std::vector<double> const& local_x,
                         std::vector<double> const& local_x_prev,
                         std::vector<double>& /*local_M_data*/,
                         std::vector<double>& /*local_K_data*/,
                         std::vector<double>& local_rhs_data,
                         std::vector<double>& local_Jac_data)
{
    auto& medium =
        *this->process_data_.media_map.getMedium(this->element_.getID());

    LocalMatrices loc_mat;
    loc_mat.setZero();
    LocalMatrices loc_mat_current_ip;
    loc_mat_current_ip.setZero();  // only to set the right matrix sizes

    typename ConstitutiveTraits::ConstitutiveSetting constitutive_setting;

    for (unsigned ip = 0; ip < this->integration_method_.getNumberOfPoints();
         ++ip)
    {
        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->element_.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, ip_data_[ip].N_u))};

        assembleWithJacobianSingleIP(
            t, dt, x_position,      //
            local_x, local_x_prev,  //
            ip_data_[ip], constitutive_setting,
            medium,              //
            loc_mat_current_ip,  //
            this->current_states_[ip], this->prev_states_[ip],
            this->material_states_[ip], this->output_data_[ip]);
        loc_mat += loc_mat_current_ip;
    }

    massLumping(loc_mat);

    addToLocalMatrixData(dt, local_x, local_x_prev, loc_mat, local_rhs_data,
                         local_Jac_data);
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    massLumping(typename ThermoRichardsMechanicsLocalAssembler<
                ShapeFunctionDisplacement, ShapeFunction, DisplacementDim,
                ConstitutiveTraits>::LocalMatrices& loc_mat) const
{
    if (this->process_data_.apply_mass_lumping)
    {
        loc_mat.storage_p_a_p =
            loc_mat.storage_p_a_p.colwise().sum().eval().asDiagonal();
        loc_mat.storage_p_a_S =
            loc_mat.storage_p_a_S.colwise().sum().eval().asDiagonal();
        loc_mat.storage_p_a_S_Jpp =
            loc_mat.storage_p_a_S_Jpp.colwise().sum().eval().asDiagonal();
    }
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    addToLocalMatrixData(
        double const dt,
        std::vector<double> const& local_x,
        std::vector<double> const& local_x_prev,
        typename ThermoRichardsMechanicsLocalAssembler<
            ShapeFunctionDisplacement, ShapeFunction, DisplacementDim,
            ConstitutiveTraits>::LocalMatrices const& loc_mat,
        std::vector<double>& local_rhs_data,
        std::vector<double>& local_Jac_data) const
{
    constexpr auto local_matrix_dim =
        displacement_size + pressure_size + temperature_size;
    assert(local_x.size() == local_matrix_dim);

    auto local_Jac = MathLib::createZeroedMatrix<
        typename ShapeMatricesTypeDisplacement::template MatrixType<
            local_matrix_dim, local_matrix_dim>>(
        local_Jac_data, local_matrix_dim, local_matrix_dim);

    auto local_rhs =
        MathLib::createZeroedVector<typename ShapeMatricesTypeDisplacement::
                                        template VectorType<local_matrix_dim>>(
            local_rhs_data, local_matrix_dim);

    local_Jac.noalias() = loc_mat.Jac;
    local_rhs.noalias() = -loc_mat.res;

    //
    // -- Jacobian
    //
    block_TT(local_Jac).noalias() += loc_mat.M_TT / dt + loc_mat.K_TT;
    block_Tp(local_Jac).noalias() +=
        loc_mat.K_Tp + loc_mat.dK_TT_dp + loc_mat.M_Tp / dt;

    block_pT(local_Jac).noalias() += loc_mat.M_pT / dt + loc_mat.K_pT;
    block_pp(local_Jac).noalias() +=
        loc_mat.K_pp + loc_mat.storage_p_a_p / dt + loc_mat.storage_p_a_S_Jpp;
    block_pu(local_Jac).noalias() = loc_mat.M_pu / dt;

    //
    // -- Residual
    //
    auto const [T, p_L, u] = localDOF(local_x);
    auto const [T_prev, p_L_prev, u_prev] = localDOF(local_x_prev);

    block_T(local_rhs).noalias() -= loc_mat.M_TT * (T - T_prev) / dt +
                                    loc_mat.K_TT * T + loc_mat.K_Tp * p_L +
                                    loc_mat.M_Tp * (p_L - p_L_prev) / dt;
    block_p(local_rhs).noalias() -=
        loc_mat.K_pp * p_L + loc_mat.K_pT * T +
        (loc_mat.storage_p_a_p + loc_mat.storage_p_a_S) * (p_L - p_L_prev) /
            dt +
        loc_mat.M_pu * (u - u_prev) / dt + loc_mat.M_pT * (T - T_prev) / dt;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    assembleWithJacobianSingleIP(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& x_position,
        std::vector<double> const& local_x,
        std::vector<double> const& local_x_prev,
        typename ThermoRichardsMechanicsLocalAssembler<
            ShapeFunctionDisplacement, ShapeFunction, DisplacementDim,
            ConstitutiveTraits>::IpData const& ip_data,
        typename ConstitutiveTraits::ConstitutiveSetting& CS,
        MaterialPropertyLib::Medium& medium,
        typename ThermoRichardsMechanicsLocalAssembler<
            ShapeFunctionDisplacement, ShapeFunction, DisplacementDim,
            ConstitutiveTraits>::LocalMatrices& out,
        typename ConstitutiveTraits::StatefulData& current_state,
        typename ConstitutiveTraits::StatefulDataPrev const& prev_state,
        MaterialStateData<DisplacementDim>& mat_state,
        typename ConstitutiveTraits::OutputData& output_data) const
{
    auto const& N_u = ip_data.N_u;
    auto const& dNdx_u = ip_data.dNdx_u;

    // N and dNdx are used for both p and T variables
    auto const& N = ip_data.N_p;
    auto const& dNdx = ip_data.dNdx_p;

    auto const B =
        LinearBMatrix::computeBMatrix<DisplacementDim,
                                      ShapeFunctionDisplacement::NPOINTS,
                                      typename BMatricesType::BMatrixType>(
            dNdx_u, N_u, (*x_position.getCoordinates())[0],
            this->is_axially_symmetric_);

    auto const [T, p_L, u] = localDOF(local_x);
    auto const [T_prev, p_L_prev, u_prev] = localDOF(local_x_prev);

    GlobalDimVectorType const grad_T_ip = dNdx * T;

    typename ConstitutiveTraits::ConstitutiveModels models(
        this->process_data_, this->solid_material_);
    typename ConstitutiveTraits::ConstitutiveTempData tmp;
    typename ConstitutiveTraits::ConstitutiveData CD;

    {
        double const T_ip = N * T;
        double const T_prev_ip = N * T_prev;

        double const p_cap_ip = -N * p_L;
        double const p_cap_prev_ip = -N * p_L_prev;
        GlobalDimVectorType const grad_p_cap_ip = -dNdx * p_L;

        KelvinVectorType eps = B * u;

        CS.eval(models, t, dt, x_position,                 //
                medium,                                    //
                {T_ip, T_prev_ip, grad_T_ip},              //
                {p_cap_ip, p_cap_prev_ip, grad_p_cap_ip},  //
                eps, current_state, prev_state, mat_state, tmp, output_data,
                CD);
    }

    using NodalMatrix = typename ShapeMatricesType::NodalMatrixType;
    NodalMatrix const NTN = N.transpose() * N;
    NodalMatrix const dNTdN = dNdx.transpose() * dNdx;

    // TODO is identity2.transpose() * B something like divergence?
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;
    typename ShapeMatricesTypeDisplacement::template MatrixType<
        displacement_size, pressure_size>
        BTI2N = B.transpose() * identity2 * N;

    /*
     * Conventions:
     *
     * * use positive signs exclusively, any negative signs must be included in
     *   the coefficients coming from the constitutive setting
     * * the used coefficients from the constitutive setting are named after the
     *   terms they appear in, e.g. K_TT_X_dNTdN means it is a scalar (X) that
     *   will be multiplied by dNdx.transpose() * dNdx. Placefolders for the
     *   coefficients are:
     *     * X -> scalar
     *     * V -> vector
     *     * K -> Kelvin vector
     * * the Laplace terms have a special name, e.g., K_TT_Laplace
     * * there shall be only one contribution to each of the LocalMatrices,
     *   assigned with = assignment; this point might be relaxed in the future
     * * this method will overwrite the data in the passed LocalMatrices& out
     *   argument, not add to it
     */

    // residual, order T, p, u
    block_p(out.res).noalias() =
        dNdx.transpose() * std::get<EqPData<DisplacementDim>>(CD).rhs_p_dNT_V;
    block_u(out.res).noalias() =
        B.transpose() *
            ProcessLib::Graph::get<TotalStressData<DisplacementDim>>(
                CD, current_state)
                .sigma_total -
        static_cast<int>(this->process_data_.apply_body_force_for_deformation) *
            N_u_op(N_u).transpose() *
            std::get<GravityData<DisplacementDim>>(CD).volumetric_body_force;

    // Storage matrices
    out.storage_p_a_p.noalias() =
        std::get<EqPData<DisplacementDim>>(CD).storage_p_a_p_X_NTN * NTN;
    out.storage_p_a_S.noalias() =
        std::get<TRMStorageData>(CD).storage_p_a_S_X_NTN * NTN;
    out.storage_p_a_S_Jpp.noalias() =
        std::get<TRMStorageData>(CD).storage_p_a_S_Jpp_X_NTN * NTN;

    // M matrices, order T, p, u
    out.M_TT.noalias() =
        std::get<EqTData<DisplacementDim>>(CD).M_TT_X_NTN * NTN;
    out.M_Tp.noalias() =
        std::get<TRMVaporDiffusionData<DisplacementDim>>(CD).M_Tp_X_NTN * NTN;

    out.M_pT.noalias() =
        std::get<EqPData<DisplacementDim>>(CD).M_pT_X_NTN * NTN;
    out.M_pu.noalias() =
        std::get<EqPData<DisplacementDim>>(CD).M_pu_X_BTI2N * BTI2N.transpose();

    // K matrices, order T, p, u
    out.K_TT.noalias() =
        dNdx.transpose() *
            std::get<TRMHeatStorageAndFluxData<DisplacementDim>>(CD)
                .K_TT_Laplace *
            dNdx +
        N.transpose() *
            (std::get<EqTData<DisplacementDim>>(CD).K_TT_NT_V_dN.transpose() *
             dNdx) +
        std::get<TRMVaporDiffusionData<DisplacementDim>>(CD).K_TT_X_dNTdN *
            dNTdN;

    out.dK_TT_dp.noalias() =
        N.transpose() *
            (std::get<TRMHeatStorageAndFluxData<DisplacementDim>>(CD)
                 .K_Tp_NT_V_dN.transpose() *
             dNdx) +
        std::get<TRMHeatStorageAndFluxData<DisplacementDim>>(CD).K_Tp_X_NTN *
            NTN;
    out.K_Tp.noalias() =
        dNdx.transpose() *
            std::get<ThermoOsmosisData<DisplacementDim>>(CD).K_Tp_Laplace *
            dNdx +
        std::get<TRMVaporDiffusionData<DisplacementDim>>(CD).K_Tp_X_dNTdN *
            dNTdN;

    out.K_pp.noalias() =
        dNdx.transpose() * std::get<EqPData<DisplacementDim>>(CD).K_pp_Laplace *
            dNdx +
        std::get<TRMVaporDiffusionData<DisplacementDim>>(CD).K_pp_X_dNTdN *
            dNTdN;
    out.K_pT.noalias() =
        dNdx.transpose() *
        std::get<ThermoOsmosisData<DisplacementDim>>(CD).K_pT_Laplace * dNdx;

    // direct Jacobian contributions, order T, p, u
    block_pT(out.Jac).noalias() =
        std::get<TRMVaporDiffusionData<DisplacementDim>>(CD).J_pT_X_dNTdN *
        dNTdN;
    block_pp(out.Jac).noalias() =
        std::get<TRMStorageData>(CD).J_pp_X_NTN * NTN +
        std::get<EqPData<DisplacementDim>>(CD).J_pp_X_BTI2NT_u_dot_N *
            BTI2N.transpose() * (u - u_prev) / dt *
            N  // TODO something with volumetric strain rate?
        + dNdx.transpose() *
              std::get<EqPData<DisplacementDim>>(CD).J_pp_dNT_V_N * N;

    block_uT(out.Jac).noalias() =
        B.transpose() *
        std::get<SolidMechanicsDataStateless<DisplacementDim>>(CD).J_uT_BT_K_N *
        N;
    block_up(out.Jac).noalias() =
        B.transpose() *
            std::get<SolidMechanicsDataStateless<DisplacementDim>>(CD)
                .J_up_BT_K_N *
            N +
        N_u_op(N_u).transpose() *
            std::get<GravityData<DisplacementDim>>(CD).J_up_HT_V_N * N;
    block_uu(out.Jac).noalias() =
        B.transpose() *
        std::get<SolidMechanicsDataStateless<DisplacementDim>>(CD)
            .stiffness_tensor *
        B;

    out *= ip_data.integration_weight;
}

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
void ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                           ShapeFunction, DisplacementDim,
                                           ConstitutiveTraits>::
    computeSecondaryVariableConcrete(double const t, double const dt,
                                     Eigen::VectorXd const& local_x,
                                     Eigen::VectorXd const& local_x_prev)
{
    auto const T = block_T(local_x);
    auto const p_L = block_p(local_x);
    auto const u = block_u(local_x);

    auto const T_prev = block_T(local_x_prev);
    auto const p_L_prev = block_p(local_x_prev);

    auto const e_id = this->element_.getID();
    auto const& process_data = this->process_data_;
    auto& medium = *process_data.media_map.getMedium(e_id);

    unsigned const n_integration_points =
        this->integration_method_.getNumberOfPoints();

    double saturation_avg = 0;
    double porosity_avg = 0;
    double liquid_density_avg = 0;
    double viscosity_avg = 0;

    using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    KV sigma_avg = KV::Zero();

    typename ConstitutiveTraits::ConstitutiveSetting constitutive_setting;

    typename ConstitutiveTraits::ConstitutiveModels models(
        process_data, this->solid_material_);
    typename ConstitutiveTraits::ConstitutiveTempData tmp;
    typename ConstitutiveTraits::ConstitutiveData CD;

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto& current_state = this->current_states_[ip];
        auto& output_data = this->output_data_[ip];

        auto const& ip_data = ip_data_[ip];

        // N is used for both p and T variables
        auto const& N = ip_data.N_p;
        auto const& N_u = ip_data.N_u;
        auto const& dNdx_u = ip_data.dNdx_u;
        auto const& dNdx = ip_data.dNdx_p;

        ParameterLib::SpatialPosition const x_position{
            std::nullopt, this->element_.getID(), ip,
            MathLib::Point3d(
                NumLib::interpolateCoordinates<ShapeFunctionDisplacement,
                                               ShapeMatricesTypeDisplacement>(
                    this->element_, N_u))};
        auto const x_coord =
            NumLib::interpolateXCoordinate<ShapeFunctionDisplacement,
                                           ShapeMatricesTypeDisplacement>(
                this->element_, N_u);
        auto const B =
            LinearBMatrix::computeBMatrix<DisplacementDim,
                                          ShapeFunctionDisplacement::NPOINTS,
                                          typename BMatricesType::BMatrixType>(
                dNdx_u, N_u, x_coord, this->is_axially_symmetric_);

        double const T_ip = N * T;
        double const T_prev_ip = N * T_prev;
        GlobalDimVectorType const grad_T_ip = dNdx * T;

        double const p_cap_ip = -N * p_L;
        double const p_cap_prev_ip = -N * p_L_prev;
        GlobalDimVectorType const grad_p_cap_ip = -dNdx * p_L;

        KelvinVectorType eps = B * u;

        constitutive_setting.eval(models,                                    //
                                  t, dt, x_position,                         //
                                  medium,                                    //
                                  {T_ip, T_prev_ip, grad_T_ip},              //
                                  {p_cap_ip, p_cap_prev_ip, grad_p_cap_ip},  //
                                  eps, current_state, this->prev_states_[ip],
                                  this->material_states_[ip], tmp, output_data,
                                  CD);

        saturation_avg += std::get<SaturationData>(current_state).S_L;
        porosity_avg += std::get<PorosityData>(current_state).phi;

        liquid_density_avg += std::get<LiquidDensityData>(output_data).rho_LR;
        viscosity_avg += std::get<LiquidViscosityData>(output_data).viscosity;
        sigma_avg += ConstitutiveTraits::ConstitutiveSetting::statefulStress(
            current_state);
    }
    saturation_avg /= n_integration_points;
    porosity_avg /= n_integration_points;
    viscosity_avg /= n_integration_points;
    liquid_density_avg /= n_integration_points;
    sigma_avg /= n_integration_points;

    (*process_data.element_saturation)[e_id] = saturation_avg;
    (*process_data.element_porosity)[e_id] = porosity_avg;
    (*process_data.element_liquid_density)[e_id] = liquid_density_avg;
    (*process_data.element_viscosity)[e_id] = viscosity_avg;

    Eigen::Map<KV>(
        &(*process_data.element_stresses)[e_id * KV::RowsAtCompileTime]) =
        MathLib::KelvinVector::kelvinVectorToSymmetricTensor(sigma_avg);

    NumLib::interpolateToHigherOrderNodes<
        ShapeFunction, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_, p_L,
                         *process_data.pressure_interpolated);
    NumLib::interpolateToHigherOrderNodes<
        ShapeFunction, typename ShapeFunctionDisplacement::MeshElement,
        DisplacementDim>(this->element_, this->is_axially_symmetric_, T,
                         *process_data.temperature_interpolated);
}
}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
