/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:27 PM
 */

#pragma once

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Property.h"

namespace ParameterLib
{
struct CoordinateSystem;
template <typename T>
struct Parameter;
}  // namespace ParameterLib

namespace MaterialPropertyLib
{
class Medium;

enum class MeanType
{
    ARITHMETIC_LINEAR,
    ARITHMETIC_SQUAREROOT,
    GEOMETRIC
};

/**
 * \brief Saturation dependent thermal conductivity model for soil.
 *
 *  The arithmetic_squareroot model is proposed by Somerton, W.~H. et al.
 *  \cite somerton1974high, which takes the form of
 *  \f[ \lambda = \lambda_{\text{dry}} + \sqrt{S}(\lambda_{\text{wet}}-
 *  \lambda_{\text{dry}}), \f]. The arithmetic_linear model is linear in
 *  the saturation and has the form \f[ \lambda = \lambda_{\text{dry}} +
 *     S(\lambda_{\text{wet}}-\lambda_{\text{dry}}),
 *  \f]
 *  The geometric model uses the weighted geometric mean
 *  \f[
 *   \lambda = \lambda_{\text{dry}}^{1-S} * \lambda_{\text{wet}}^S,
 *  \f]
 *  where \f$\lambda_{\text{dry}}\f$ is the thermal conductivity of soil at the
 *  dry state, \f$\lambda_{\text{wet}}\f$ is the thermal conductivity of soil at
 *  the fully water saturated state, and \f$S\f$ is the water saturation.
 */
template <MeanType MeantType, int GlobalDimension>
class SaturationWeightedThermalConductivity final : public Property
{
public:
    SaturationWeightedThermalConductivity(
        std::string name,
        ParameterLib::Parameter<double> const& dry_thermal_conductivity,
        ParameterLib::Parameter<double> const& wet_thermal_conductivity);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationWeightedThermalConductivity' is "
                "implemented on the 'media' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t,
                            double const dt) const override;

private:
    /// Thermal conductivity of soil at the dry state.
    ParameterLib::Parameter<double> const& dry_thermal_conductivity_;
    /// Thermal conductivity of soil at the fully water saturated state.
    ParameterLib::Parameter<double> const& wet_thermal_conductivity_;
};

extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_LINEAR,
    1>;
extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_SQUAREROOT,
    1>;
extern template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC,
                                                            1>;
extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_LINEAR,
    2>;
extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_SQUAREROOT,
    2>;
extern template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC,
                                                            2>;
extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_LINEAR,
    3>;
extern template class SaturationWeightedThermalConductivity<
    MeanType::ARITHMETIC_SQUAREROOT,
    3>;
extern template class SaturationWeightedThermalConductivity<MeanType::GEOMETRIC,
                                                            3>;

}  // namespace MaterialPropertyLib
