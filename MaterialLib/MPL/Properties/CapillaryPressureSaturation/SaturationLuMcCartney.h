// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <boost/math/differentiation/autodiff.hpp>
#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

#include "MaterialLib/MPL/Properties/GibbsFreeEnergy/DimensionLessGibbsFreeEnergyRegion1.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \brief The Lu - Mc Cartney water retention behavior.
 * \details This property must be a medium property, it
 * computes the saturation of the wetting phase as function
 * of capillary pressure.
 *
 * Original source: Lu, Y., McCartney, J.S. Temperature effects on adsorption
and capillarity water retention mechanisms in constrained unsaturated soils.
Acta Geotech. 19, 6467–6482 (2024).
 * https://doi.org/10.1007/s11440-024-02341-9
 * and Guo, G., Zheng, L., Lu, Y., Behbehani, F., & McCartney, J. (2026).
Coupled THM modeling of bentonite heating and hydration in tank TemperatureTypes
with a new temperature-dependent water retention model. Computers and
Geotechnics, 190, 107736.
 * https://doi.org/10.1016/j.compgeo.2025.107736
Get rights and content

 */
class SaturationLuMcCartney final : public Property
{
public:
    SaturationLuMcCartney(std::string name, std::string const& material);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationLuMcCartney' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/,
                           double const /*dt*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/,
                            double const /*dt*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1, Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/,
                             double const /*dt*/) const override;

private:
    template <typename PressureType, typename TemperatureType>
    boost::math::differentiation::promote<PressureType, TemperatureType>
    water_content(PressureType p_cap, TemperatureType T,
                  VariableArray const& variable_array) const;
    template <typename PressureType, typename TemperatureType>
    boost::math::differentiation::promote<PressureType, TemperatureType>
    adsorptive_water_content(PressureType p_cap, TemperatureType T,
                             VariableArray const& variable_array) const;
    template <typename PressureType, typename TemperatureType>
    boost::math::differentiation::promote<PressureType, TemperatureType>
    capillary_water_content(PressureType p_cap, TemperatureType T,
                            VariableArray const& variable_array) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> theta_a_max(
        TemperatureType T, VariableArray const& variable_array) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> psi_max(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> c(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> CEC(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> psi_c(
        TemperatureType T, VariableArray const& variable_array) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> A_H(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> chi(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> delta_h(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> alpha(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> epsilon_w(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> zeta_s(
        TemperatureType T) const;
    template <typename TemperatureType>
    boost::math::differentiation::promote<TemperatureType> rho_w(
        TemperatureType T) const;

    std::string const material_;

    const MaterialLib::Fluid::DimensionLessGibbsFreeEnergyRegion1
        gibbs_free_energy_;
    // constants
    double const R = 8.31432;         ///< universal gas constant in J/(mol K)
    double const kB = 1.3806505e-23;  ///< Boltzmann constant in J/K
    double const hp = 6.626068e-34;   ///< Planck constant in J s
    static constexpr double ref_T_ = 1386;  ///< reference temperature in K.
    static constexpr double ref_p_ = 1.653e7;

    /// <---- Lu and McCartney 2024 parameters ---->
    double const nu_w = 1.8e-5;  ///< molar volume of water in m^3/mol
    double const A = 0.92;  ///< fitting parameter for cation exchange capacity
    double const B = 0.9;   ///< fitting parameter for cation exchange capacity
    double const T1 = 373.15;  ///< parameter for temperature effect of CEC in K
    double const T2 =
        1273.15;  ///< parameter for temperature effect of CEC in K
    double const nu_e =
        2.45e9;                ///< main absorption frequency of UV light in Hz
    double const Tr = 293.15;  ///< reference temperature in K
    double const C1 = -0.00151;  ///< parameter for temperature correction
    double const n_w = 1.33;     ///< diffraction index

    /// material parameters for adsorptive water content
    double E1_minus_EL;  ///< energy difference between heat of adorption
                         ///< of first and second layer in J/mol
    double nu_mr;        ///< gravimetric moisture content when first layer is
                         ///< fully saturated
    double CEC_max;      ///< maximum cation exchange capacity in meq/g
    double bw;
    double M;  ///< exponent for adsorptive water content

    /// material parameters for capillary water content
    double SSA;        ///< specific surface area in m^2/g
    double epsilon_s;  ///<  dielectic constant of soil particles
    double n_s;        ///< refractive index of soil particles
    double eta_alpha;  ///< soil related coefficient
    double
        alpha_0_inv;  ///< air entry suction at reference temperature Tr in Pa
    double chi_r;
    double N;  ///< exponent for capillary water content
};
}  // namespace MaterialPropertyLib
