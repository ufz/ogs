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

#include <Eigen/Eigen>
#include <memory>
#include <utility>

#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace PhaseField
{
enum class PhaseFieldModel
{
    AT1,
    AT2,
    COHESIVE
};

enum class SofteningCurve
{
    Linear,
    Exponential
};

enum class EnergySplitModel
{
    Isotropic,
    VolDev,
    EffectiveStress
};

class DegradationDerivative
{
public:
    virtual ~DegradationDerivative() = default;
    virtual double degradation(double const d_ip,
                               double const k,
                               double const ls) = 0; /* degradation_df0 */
    virtual double degradation_df1(double const d_ip, double const ls) = 0;
    virtual double degradation_df2(double const d_ip, double const ls) = 0;
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
    double degradation_df1(double const d_ip, double const /* ls */) override
    {
        return 2. * d_ip;
    };
    double degradation_df2(double const /* d_ip */,
                           double const /* ls */) override
    {
        return 2.;
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
    double degradation_df1(double const d_ip, double const ls) override
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
                return (n * std::pow(d_ip, n - 1.) * a1 -
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
                return (n * std::pow(d_ip, n - 1.) * a1 -
                        std::pow(d_ip, n) * a2) /
                       (a1 * a1);
            }
        }
    };
    double degradation_df2(double const d_ip, double const ls) override
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
                return (2. * a2 * a2 * std::pow(d_ip, n) -
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
                return (2. * a2 * a2 * std::pow(d_ip, n) -
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
struct PhaseFieldProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    ParameterLib::Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    bool hydro_crack = false;
    bool crack_pressure = false;
    double irreversible_threshold;
    PhaseFieldModel phasefield_model;
    EnergySplitModel energy_split_model;
    SofteningCurve softening_curve;
    double characteristic_length;
    std::unique_ptr<DegradationDerivative> degradation_derivative;

    double const unity_pressure = 1.0;
    double pressure = 0.0;
    double pressure_old = 0.0;
    double pressure_error = 0.0;
    double injected_volume = 0.0;
    double crack_volume = 0.0;
    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
};

}  // namespace PhaseField
}  // namespace ProcessLib
