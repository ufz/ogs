/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
/**
 * \brief A saturation dependent thermal conductivity model for soil.
 *
 *  This model is proposed by Somerton, W.~H. et al. \cite somerton1974high,
 *  which takes the form of
 *  \f[
 *   \lambda = \lambda_{\text{dry}} +
 *     \sqrt{S}(\lambda_{\text{wet}}-\lambda_{\text{dry}}),
 *  \f]
 * where \f$\lambda_{\text{dry}}\f$ is the thermal conductivity of soil at the
 * dry state, \f$\lambda_{\text{wet}}\f$ is the thermal conductivity of soil at
 * the fully water saturated state, and \f$S\f$ is the water saturation.
 */
template <int GlobalDimension>
class SoilThermalConductivitySomerton final : public Property
{
public:
    SoilThermalConductivitySomerton(
        std::string name,
        ParameterLib::Parameter<double> const& dry_thermal_conductivity,
        ParameterLib::Parameter<double> const& wet_thermal_conductivity,
        ParameterLib::CoordinateSystem const* const local_coordinate_system);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SoilThermalConductivitySomerton' is "
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

    ParameterLib::CoordinateSystem const* const local_coordinate_system_;
};

extern template class SoilThermalConductivitySomerton<2>;
extern template class SoilThermalConductivitySomerton<3>;

}  // namespace MaterialPropertyLib
