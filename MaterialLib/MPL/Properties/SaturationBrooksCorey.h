/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include <limits>
#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \brief A well known soil characteristics function
 * \details This property must be a medium property, it
 * computes the saturation of the wetting phase as function
 * of capillary pressure.
 */
class SaturationBrooksCorey final : public Property
{
private:
    Medium* _medium = nullptr;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _exponent;
    const double _entry_pressure;

public:
    SaturationBrooksCorey(const double residual_liquid_saturation,
                          const double residual_gas_saturation,
                          const double exponent,
                          const double entry_pressure);

    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (!std::holds_alternative<Medium*>(scale_pointer))
        {
            OGS_FATAL(
                "The property 'SaturationBrooksCorey' is implemented on the "
                "'media' scale only.");
        }
        _medium = std::get<Medium*>(scale_pointer);
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& /*pos*/,
                           double const /*t*/) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& /*pos*/,
                            double const /*t*/) const override;
    PropertyDataType d2Value(VariableArray const& variable_array,
                             Variable const variable1,
                             Variable const variable2,
                             ParameterLib::SpatialPosition const& /*pos*/,
                             double const /*t*/) const override;
};

inline std::unique_ptr<SaturationBrooksCorey> createSaturationBrooksCorey(
    BaseLib::ConfigTree const& config)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey}
    // config.checkConfigParameter("type", "SaturationBrooksCorey");
    DBUG("Create SaturationBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const exponent =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    auto const entry_pressure =
        //! \ogs_file_param{prj__media__medium__properties__property__SaturationBrooksCorey__entry_pressure}
        config.getConfigParameter<double>("entry_pressure");

    return std::make_unique<SaturationBrooksCorey>(residual_liquid_saturation,
                                                   residual_gas_saturation,
                                                   exponent, entry_pressure);
}

}  // namespace MaterialPropertyLib
