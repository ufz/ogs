// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <numbers>
#include <optional>

#include "BaseLib/Logging.h"
#include "LinearElasticIsotropic.h"
#include "LinearElasticIsotropicPhaseField.h"
#include "LinearElasticOrthotropic.h"
#include "LinearElasticOrthotropicPhaseField.h"

namespace MaterialLib
{
namespace Solids
{
namespace Phasefield
{
enum class PhaseFieldModel
{
    AT1,
    AT2,
    COHESIVE
};

template <int DisplacementDim>
PhaseFieldModel convertStringToPhaseFieldModel(
    std::string const& phasefield_model)
{
    if (phasefield_model == "AT1")
    {
        return PhaseFieldModel::AT1;
    }
    if (phasefield_model == "AT2")
    {
        return PhaseFieldModel::AT2;
    }
    if (phasefield_model == "COHESIVE")
    {
        return PhaseFieldModel::COHESIVE;
    }
    OGS_FATAL(
        "phasefield_model must be 'AT1', 'AT2' or 'COHESIVE' but '{}' "
        "was given",
        phasefield_model.c_str());
};

enum class SofteningCurve
{
    Linear,
    Exponential
};

template <int DisplacementDim>
SofteningCurve convertStringToSofteningCurve(
    std::optional<std::string> const& softening_curve)
{
    if (softening_curve)
    {
        if (*softening_curve == "Linear")
        {
            return SofteningCurve::Linear;
        }
        if (*softening_curve == "Exponential")
        {
            return SofteningCurve::Exponential;
        }
        OGS_FATAL(
            "softening_curve must be 'Linear' or 'Exponential' but '{}' "
            "was given",
            softening_curve->c_str());
    }
    return SofteningCurve::Linear;  // default
};

enum class EnergySplitModel
{
    Isotropic,
    VolDev,
    Spectral,
    EffectiveStress,
    OrthoVolDev,
    OrthoMasonry
};

template <int DisplacementDim>
EnergySplitModel convertStringToEnergySplitModel(
    std::string const& energy_split_model)
{
    if (energy_split_model == "Isotropic")
    {
        return EnergySplitModel::Isotropic;
    }
    if (energy_split_model == "VolumetricDeviatoric")
    {
        return EnergySplitModel::VolDev;
    }
    if (energy_split_model == "Spectral")
    {
        return EnergySplitModel::Spectral;
    }
    if (energy_split_model == "EffectiveStress")
    {
        return EnergySplitModel::EffectiveStress;
    }
    if (energy_split_model == "OrthoVolDev")
    {
        return EnergySplitModel::OrthoVolDev;
    }
    if (energy_split_model == "OrthoMasonry")
    {
        return EnergySplitModel::OrthoMasonry;
    }
    OGS_FATAL(
        "energy_split_model must be 'Isotropic', 'VolumetricDeviatoric', "
        "'EffectiveStress', 'OrthoVolDev' or 'OrthoMasonry' but '{}' was given",
        energy_split_model.c_str());
};

class DegradationDerivative
{
public:
    virtual ~DegradationDerivative() = default;
    virtual double degradation(double const d_ip,
                               double const k,
                               double const ls) = 0; /* degradation_df0 */
    virtual double degradationDf1(double const d_ip, double const k,
                                  double const ls) = 0;
    virtual double degradationDf2(double const d_ip, double const k,
                                  double const ls) = 0;
};

class AT_DegradationDerivative : public DegradationDerivative
{
public:
    double degradation(double const d_ip,
                       double const k,
                       double const /* ls */) override
    {
        return d_ip * d_ip * (1. - k) + k;
    };
    double degradationDf1(double const d_ip, double const k,
                          double const /* ls */) override
    {
        return 2. * (1. - k) * d_ip;
    };
    double degradationDf2(double const /* d_ip */, double const k,
                          double const /* ls */) override
    {
        return 2. * (1. - k);
    };
};

class COHESIVE_DegradationDerivative : public DegradationDerivative
{
private:
    double const lch;
    SofteningCurve softening_curve;

public:
    COHESIVE_DegradationDerivative(double const Parameter_lch,
                                   SofteningCurve Parameter_softening_curve)
        : lch(Parameter_lch), softening_curve(Parameter_softening_curve){};
    double degradation(double const d_ip,
                       double const k,
                       double const ls) override
    {
        double const m1 = 4.0 * lch / acos(-1.) / ls;
        switch (softening_curve)
        {
            case SofteningCurve::Exponential:
            {
                double const m2 = std::pow(2., 5. / 3.) - 3.;
                double const n = 2.5;
                return std::pow(d_ip, n) /
                           (std::pow(d_ip, n) +
                            m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip))) *
                           (1.0 - k) +
                       k;
            }
            default:
            {
                double const m2 = -0.5;
                double const n = 2.;
                return std::pow(d_ip, n) /
                           (std::pow(d_ip, n) +
                            m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip))) *
                           (1. - k) +
                       k;
            }
        }
    };
    double degradationDf1(double const d_ip, double const k,
                          double const ls) override
    {
        double const m1 = 4.0 * lch / acos(-1.) / ls;
        switch (softening_curve)
        {
            case SofteningCurve::Exponential:
            {
                double const m2 = std::pow(2., 5. / 3.) - 3.;
                double const n = 2.5;
                double const a1 = std::pow(d_ip, n) +
                                  m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip));
                double const a2 = n * std::pow(d_ip, n - 1.) -
                                  2. * m1 * m2 * (1. - d_ip) - m1;
                return (1. - k) *
                       (n * std::pow(d_ip, n - 1.) * a1 -
                        std::pow(d_ip, n) * a2) /
                       (a1 * a1);
            }
            default:
            {
                double const m2 = -0.5;
                double const n = 2.;
                double const a1 = std::pow(d_ip, n) +
                                  m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip));
                double const a2 = n * std::pow(d_ip, n - 1.) -
                                  2. * m1 * m2 * (1. - d_ip) - m1;
                return (1. - k) *
                       (n * std::pow(d_ip, n - 1.) * a1 -
                        std::pow(d_ip, n) * a2) /
                       (a1 * a1);
            }
        }
    };
    double degradationDf2(double const d_ip, double const k,
                          double const ls) override
    {
        double const m1 = 4.0 * lch / acos(-1.) / ls;
        switch (softening_curve)
        {
            case SofteningCurve::Exponential:
            {
                double const m2 = std::pow(2., 5. / 3.) - 3.;
                double const n = 2.5;
                double a1 = std::pow(d_ip, n) +
                            m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip));
                double a2 = n * std::pow(d_ip, n - 1.) -
                            2. * m1 * m2 * (1. - d_ip) - m1;
                return (1. - k) *
                       (2. * a2 * a2 * std::pow(d_ip, n) -
                        a1 * std::pow(d_ip, n) *
                            (2. * m1 * m2 +
                             n * std::pow(d_ip, n - 2.) * (n - 1.)) -
                        2. * a1 * a2 * n * std::pow(d_ip, n - 1.) +
                        a1 * a1 * n * std::pow(d_ip, n - 2.) * (n - 1.)) /
                       (std::pow(a1, 3));
            }
            default:
            {
                double const m2 = -0.5;
                double const n = 2.;
                double const a1 = std::pow(d_ip, n) +
                                  m1 * (1. - d_ip) * (1. + m2 * (1. - d_ip));
                double const a2 = n * std::pow(d_ip, n - 1.) -
                                  2. * m1 * m2 * (1. - d_ip) - m1;
                return (1. - k) *
                       (2. * a2 * a2 * std::pow(d_ip, n) -
                        a1 * std::pow(d_ip, n) *
                            (2. * m1 * m2 +
                             n * std::pow(d_ip, n - 2.) * (n - 1.)) -
                        2. * a1 * a2 * n * std::pow(d_ip, n - 1.) +
                        a1 * a1 * n * std::pow(d_ip, n - 2.) * (n - 1.)) /
                       (std::pow(a1, 3));
            }
        }
    };
};

template <int DisplacementDim>
std::unique_ptr<DegradationDerivative> creatDegradationDerivative(
    PhaseFieldModel const& phasefield_model, double const lch,
    SofteningCurve const& softening_curve)
{
    switch (phasefield_model)
    {
        case PhaseFieldModel::COHESIVE:
            return std::make_unique<COHESIVE_DegradationDerivative>(
                lch, softening_curve);
        default:
            return std::make_unique<AT_DegradationDerivative>();
    }
};

template <typename T_DNDX, typename T_N, typename T_W, typename T_D,
          typename T_LOCAL_JAC, typename T_LOCAL_RHS>
void calculateCrackLocalJacobianAndResidual(T_DNDX& dNdx, T_N& N, T_W& w,
                                            T_D& d, T_LOCAL_JAC& local_Jac,
                                            T_LOCAL_RHS local_rhs,
                                            double const gc, double const ls,
                                            PhaseFieldModel& phasefield_model)
{
    switch (phasefield_model)
    {
        case PhaseFieldModel::AT1:
        {
            auto const local_Jac_AT1 =
                (gc * 0.75 * dNdx.transpose() * dNdx * ls * w).eval();
            local_Jac.noalias() += local_Jac_AT1;

            local_rhs.noalias() -=
                gc * (-0.375 * N.transpose() / ls) * w + local_Jac_AT1 * d;
            break;
        }
        case PhaseFieldModel::AT2:
        {
            auto const local_Jac_AT2 =
                (gc * (N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) *
                 w)
                    .eval();
            local_Jac.noalias() += local_Jac_AT2;

            local_rhs.noalias() -=
                local_Jac_AT2 * d - gc * (N.transpose() / ls) * w;
            break;
        }
        case PhaseFieldModel::COHESIVE:
        {
            auto const local_Jac_COHESIVE =
                (2.0 / std::numbers::pi * gc *
                 (-N.transpose() * N / ls + dNdx.transpose() * dNdx * ls) * w)
                    .eval();

            local_Jac.noalias() += local_Jac_COHESIVE;

            local_rhs.noalias() -= local_Jac_COHESIVE * d;
            break;
        }
        default:
        {
            OGS_FATAL("Invalid phase field model specified.");
        }
    }
};

template <typename T_VECTOR, typename T_MATRIX, int DisplacementDim>
void calculateStress(
    T_VECTOR& sigma, T_VECTOR& sigma_tensile, T_VECTOR& sigma_compressive,
    T_VECTOR& eps_tensile, T_MATRIX& D, T_MATRIX& C_tensile,
    T_MATRIX& C_compressive, double& strain_energy_tensile,
    double& elastic_energy, double const degradation, T_VECTOR const& eps,
    EnergySplitModel const& energy_split_model, double const t,
    ParameterLib::SpatialPosition const& x,
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material)
{
    auto linear_elastic_mp =
        static_cast<MaterialLib::Solids::LinearElasticIsotropic<
            DisplacementDim> const&>(solid_material)
            .getMaterialProperties();

    auto const lambda = linear_elastic_mp.lambda(t, x);
    auto const bulk_modulus = linear_elastic_mp.bulk_modulus(t, x);
    auto const mu = linear_elastic_mp.mu(t, x);
    switch (energy_split_model)
    {
        case EnergySplitModel::Isotropic:
        {
            std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                     elastic_energy, C_tensile, C_compressive) =
                MaterialLib::Solids::Phasefield::
                    calculateIsotropicDegradedStress<DisplacementDim>(
                        degradation, bulk_modulus, mu, eps);
            break;
        }
        case EnergySplitModel::VolDev:
        {
            std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                     elastic_energy, C_tensile, C_compressive) =
                MaterialLib::Solids::Phasefield::calculateVolDevDegradedStress<
                    DisplacementDim>(degradation, bulk_modulus, mu, eps);
            break;
        }
        case EnergySplitModel::Spectral:
        {
            std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                     elastic_energy, C_tensile, C_compressive) =
                MaterialLib::Solids::Phasefield::
                    calculateSpectralDegradedStress<DisplacementDim>(
                        degradation, lambda, mu, eps);
            break;
        }
        case EnergySplitModel::EffectiveStress:
        {
            std::tie(sigma, sigma_tensile, D, strain_energy_tensile,
                     elastic_energy, C_tensile, C_compressive) =
                MaterialLib::Solids::Phasefield::
                    calculateIsotropicDegradedStressWithRankineEnergy<
                        DisplacementDim>(degradation, bulk_modulus, mu, eps);
            break;
        }
        case EnergySplitModel::OrthoVolDev:
        {
            double temperature = 0.;
            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C_ortho =
                static_cast<MaterialLib::Solids::LinearElasticOrthotropic<
                    DisplacementDim> const&>(solid_material)
                    .getElasticTensor(t, x, temperature);

            std::tie(eps_tensile, sigma, sigma_tensile, sigma_compressive, D,
                     strain_energy_tensile, elastic_energy, C_tensile,
                     C_compressive) = MaterialLib::Solids::Phasefield::
                calculateOrthoVolDevDegradedStress<DisplacementDim>(
                    degradation, eps, C_ortho);
            break;
        }
        case EnergySplitModel::OrthoMasonry:
        {
            double temperature = 0.;
            MathLib::KelvinVector::KelvinMatrixType<DisplacementDim> C_ortho =
                static_cast<MaterialLib::Solids::LinearElasticOrthotropic<
                    DisplacementDim> const&>(solid_material)
                    .getElasticTensor(t, x, temperature);

            std::tie(eps_tensile, sigma, sigma_tensile, D,
                     strain_energy_tensile, elastic_energy, C_tensile,
                     C_compressive) = MaterialLib::Solids::Phasefield::
                calculateOrthoMasonryDegradedStress<DisplacementDim>(
                    degradation, eps, C_ortho);
            break;
        }
        default:
            OGS_FATAL("Invalid energy split model specified.");
    }
};
}  // namespace Phasefield
}  // namespace Solids
}  // namespace MaterialLib
