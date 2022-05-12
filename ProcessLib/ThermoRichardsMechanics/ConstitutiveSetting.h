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

#include <limits>

#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MathLib/KelvinVector.h"
#include "ProcessLib/ThermoRichardsMechanics/ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
static constexpr double nan = std::numeric_limits<double>::quiet_NaN();

template <int DisplacementDim>
struct ConstitutiveSettingEqCommon
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using GlobalDimVectorType = Eigen::Vector<double, DisplacementDim>;
    using GlobalDimMatrixType = Eigen::Matrix<double, DisplacementDim,
                                              DisplacementDim, Eigen::RowMajor>;

    void eval(double const t, double const dt,
              ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium, double const T,
              double const T_dot, GlobalDimVectorType const grad_T,
              double const p_cap, double const p_cap_dot,
              GlobalDimVectorType const grad_p_cap, KelvinVector const& eps_arg,
              KelvinVector const& eps_prev_arg);
};

template <int DisplacementDim>
struct ConstitutiveSettingEqU
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using GlobalDimVectorType = Eigen::Vector<double, DisplacementDim>;
    using GlobalDimMatrixType = Eigen::Matrix<double, DisplacementDim,
                                              DisplacementDim, Eigen::RowMajor>;

    explicit ConstitutiveSettingEqU(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data)
        : solid_material_(solid_material),
          material_state_variables_(
              solid_material.createMaterialStateVariables()),
          process_data_(process_data)
    {
        // Initialize current time step values
        sigma_eff.setZero(kelvin_vector_size);
        sigma_sw.setZero(kelvin_vector_size);
        eps.setZero(kelvin_vector_size);
        eps_m.setZero(kelvin_vector_size);

        // Previous time step values are not initialized and are set later.
        eps_prev.resize(kelvin_vector_size);
        eps_m_prev.resize(kelvin_vector_size);
        sigma_eff_prev.resize(kelvin_vector_size);
    }

    void eval(double const t, double const dt,
              ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium, double const T,
              double const T_dot, GlobalDimVectorType const& grad_T,
              double const p_cap, double const p_cap_dot,
              GlobalDimVectorType const& grad_p_cap,
              KelvinVector const& eps_arg, KelvinVector const& eps_prev_arg);

    void pushBackState()
    {
        eps_prev = eps;
        eps_m_prev = eps_m;
        sigma_eff_prev = sigma_eff;
        sigma_sw_prev = sigma_sw;
        material_state_variables_->pushBackState();
    }

    KelvinMatrix computeElasticTangentStiffness(
        double const t, ParameterLib::SpatialPosition const& x_position,
        double const dt, double const temperature_prev,
        double const temperature)
    {
        namespace MPL = MaterialPropertyLib;

        MPL::VariableArray variable_array;
        MPL::VariableArray variable_array_prev;

        auto const null_state = solid_material_.createMaterialStateVariables();

        using KV = MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;

        variable_array[static_cast<int>(MPL::Variable::stress)].emplace<KV>(
            KV::Zero());
        variable_array[static_cast<int>(MPL::Variable::mechanical_strain)]
            .emplace<KV>(KV::Zero());
        variable_array[static_cast<int>(MPL::Variable::temperature)]
            .emplace<double>(temperature);

        variable_array_prev[static_cast<int>(MPL::Variable::stress)]
            .emplace<KV>(KV::Zero());
        variable_array_prev[static_cast<int>(MPL::Variable::mechanical_strain)]
            .emplace<KV>(KV::Zero());
        variable_array_prev[static_cast<int>(MPL::Variable::temperature)]
            .emplace<double>(temperature_prev);

        auto&& solution =
            solid_material_.integrateStress(variable_array_prev, variable_array,
                                            t, x_position, dt, *null_state);

        if (!solution)
        {
            OGS_FATAL("Computation of elastic tangent stiffness failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C =
            std::move(std::get<2>(*solution));

        return C;
    }

    KelvinMatrix updateConstitutiveRelation(
        MaterialPropertyLib::VariableArray const& variable_array,
        double const t,
        ParameterLib::SpatialPosition const& x_position,
        double const dt,
        double const temperature_prev)
    {
        MaterialPropertyLib::VariableArray variable_array_prev;
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::stress)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                sigma_eff_prev);
        variable_array_prev[static_cast<int>(MaterialPropertyLib::Variable::
                                                 mechanical_strain)]
            .emplace<MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                eps_m_prev);
        variable_array_prev[static_cast<int>(
                                MaterialPropertyLib::Variable::temperature)]
            .emplace<double>(temperature_prev);

        auto solution = solid_material_.integrateStress(
            variable_array_prev, variable_array, t, x_position, dt,
            *material_state_variables_);

        if (!solution)
        {
            OGS_FATAL("Computation of local constitutive relation failed.");
        }

        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C;
        std::tie(sigma_eff, material_state_variables_, C) =
            std::move(*solution);

        return C;
    }

    static KelvinVector KVnan()
    {
        return KelvinVector::Constant(kelvin_vector_size, nan);
    }
    static KelvinMatrix KMnan()
    {
        return KelvinMatrix::Constant(kelvin_vector_size, kelvin_vector_size,
                                      nan);
    }
    static GlobalDimVectorType DVnan()
    {
        return GlobalDimVectorType::Constant(kelvin_vector_size, nan);
    }

    KelvinVector sigma_eff = KVnan();
    KelvinVector sigma_eff_prev = KVnan();
    KelvinVector sigma_sw = KVnan();
    KelvinVector sigma_sw_prev = KVnan();
    KelvinVector eps_m = KVnan();
    KelvinVector eps_m_prev = KVnan();
    KelvinVector eps = KVnan();
    KelvinVector eps_prev = KVnan();

    KelvinVector J_up_BT_K_N = KVnan();
    GlobalDimVectorType J_up_HT_V_N = DVnan();
    double J_up_X_BTI2N = nan;
    KelvinVector J_uT_BT_K_N = KVnan();

    KelvinVector sigma_total = KVnan();
    KelvinMatrix stiffness_tensor = KMnan();

    GlobalDimVectorType volumetric_body_force = DVnan();

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariables() const
    {
        return *material_state_variables_;
    }

    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material_;
    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables_;

    ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data_;
};

template <int DisplacementDim>
struct ConstitutiveSettingEqP
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using GlobalDimVectorType = Eigen::Vector<double, DisplacementDim>;
    using GlobalDimMatrixType = Eigen::Matrix<double, DisplacementDim,
                                              DisplacementDim, Eigen::RowMajor>;

    void eval(double const t, double const dt,
              ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium, double const T,
              double const T_dot, GlobalDimVectorType const& grad_T,
              double const p_cap, double const p_cap_dot,
              GlobalDimVectorType const& grad_p_cap,
              KelvinVector const& eps_arg, KelvinVector const& eps_prev_arg);

    static GlobalDimVectorType DVnan()
    {
        return GlobalDimVectorType::Constant(kelvin_vector_size, nan);
    }
    static GlobalDimMatrixType DMnan()
    {
        return GlobalDimMatrixType::Constant(kelvin_vector_size,
                                             kelvin_vector_size, nan);
    }

    GlobalDimVectorType J_pp_dNT_V_N = DVnan();
    double J_pp_X_BTI2NT_u_dot_N = nan;
    double J_pp_X_NTN = nan;
    double J_pT_X_dNTdN = nan;

    GlobalDimMatrixType K_pp_Laplace = DMnan();
    double K_pp_X_dNTdN = nan;

    double M_pT_X_NTN = nan;
    double M_pu_X_BTI2N = nan;

    GlobalDimVectorType rhs_p_dNT_V = DVnan();

    double storage_p_a_p_X_NTN = nan;
    double storage_p_a_S_Jpp_X_NTN = nan;
    double storage_p_a_S_X_NTN = nan;
};

template <int DisplacementDim>
struct ConstitutiveSettingEqT
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using GlobalDimVectorType = Eigen::Vector<double, DisplacementDim>;
    using GlobalDimMatrixType = Eigen::Matrix<double, DisplacementDim,
                                              DisplacementDim, Eigen::RowMajor>;

    void eval(double const t, double const dt,
              ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium, double const T,
              double const T_dot, GlobalDimVectorType const& grad_T,
              double const p_cap, double const p_cap_dot,
              GlobalDimVectorType const& grad_p_cap,
              KelvinVector const& eps_arg, KelvinVector const& eps_prev_arg);

    static GlobalDimVectorType DVnan()
    {
        return GlobalDimVectorType::Constant(kelvin_vector_size, nan);
    }
    static GlobalDimMatrixType DMnan()
    {
        return GlobalDimMatrixType::Constant(kelvin_vector_size,
                                             kelvin_vector_size, nan);
    }

    GlobalDimVectorType K_Tp_NT_V_dN = DVnan();
    double K_Tp_X_dNTdN = nan;
    double K_Tp_X_NTN = nan;
    GlobalDimMatrixType K_TT_Laplace = DMnan();
    GlobalDimVectorType K_TT_NT_V_dN = DVnan();
    double K_TT_X_dNTdN = nan;

    double M_Tp_X_NTN = nan;
    double M_TT_X_NTN = nan;
};

template <int DisplacementDim>
struct ConstitutiveSetting
{
    static constexpr int kelvin_vector_size =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<kelvin_vector_size>;

    using KelvinVector =
        MathLib::KelvinVector::KelvinVectorType<DisplacementDim>;
    using KelvinMatrix =
        MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>;

    using GlobalDimVectorType = Eigen::Vector<double, DisplacementDim>;
    using GlobalDimMatrixType = Eigen::Matrix<double, DisplacementDim,
                                              DisplacementDim, Eigen::RowMajor>;

    explicit ConstitutiveSetting(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material,
        ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data)
        : eqU(solid_material, process_data), process_data_(process_data)
    {
    }

    void eval(double const t, double const dt,
              ParameterLib::SpatialPosition const& x_position,
              MaterialPropertyLib::Medium& medium, double const T,
              double const T_dot, GlobalDimVectorType const& grad_T,
              double const p_cap, double const p_cap_dot,
              GlobalDimVectorType const& grad_p_cap,
              KelvinVector const& eps_arg, KelvinVector const& eps_prev_arg);

    void pushBackState()
    {
        saturation_prev = saturation;
        porosity_prev = porosity;
        transport_porosity_prev = transport_porosity;

        eqU.pushBackState();
    }

    static GlobalDimVectorType DVnan()
    {
        return GlobalDimVectorType::Constant(kelvin_vector_size, nan);
    }
    GlobalDimVectorType v_darcy = DVnan();

    double saturation = nan;
    double saturation_prev = nan;
    double porosity = nan;
    double porosity_prev = nan;
    double transport_porosity = nan;
    double transport_porosity_prev = nan;
    double liquid_density = nan;
    double viscosity = nan;
    double dry_density_solid = nan;

    ConstitutiveSettingEqU<DisplacementDim> eqU;
    ConstitutiveSettingEqP<DisplacementDim> eqP;
    ConstitutiveSettingEqT<DisplacementDim> eqT;

private:
    ConstitutiveSettingEqCommon<DisplacementDim> com_;

    ThermoRichardsMechanicsProcessData<DisplacementDim> const& process_data_;
};

extern template struct ConstitutiveSetting<2>;
extern template struct ConstitutiveSetting<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
