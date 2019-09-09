/**
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

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class RelPermBrooksCorey
 * \brief Relative permeability function proposed by Brooks&Corey
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function
 * of capillary pressure.
 */
class RelPermBrooksCorey final : public Property
{
private:
    Medium* _medium;
    const double _residual_liquid_saturation;
    const double _residual_gas_saturation;
    const double _min_relative_permeability_liquid;
    const double _min_relative_permeability_gas;
    const double _exponent;

public:
    RelPermBrooksCorey(const double /*residual_liquid_saturation*/,
                       const double /*residual_gas_saturation*/,
                       const double /*_min_relative_permeability_liquid*/,
                       const double /*_min_relative_permeability_gas*/,
                       const double /*exponent*/
    );
    /// This method assigns a pointer to the meterial object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Medium*>(scale_pointer))
        {
            _medium = std::get<Medium*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'RelPermBrooksCorey' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t) const override;
};

inline std::unique_ptr<RelPermBrooksCorey> createRelPermBrooksCorey(
    BaseLib::ConfigTree const& config)
{
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey}
    // config.checkConfigParameter("type", "RelPermBrooksCorey");
    DBUG("Create RelPermBrooksCorey medium property");

    auto const residual_liquid_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_liquid_saturation}
        config.getConfigParameter<double>("residual_liquid_saturation");
    auto const residual_gas_saturation =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__residual_gas_saturation}
        config.getConfigParameter<double>("residual_gas_saturation");
    auto const min_relative_permeability_liquid =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_liquid}
        config.getConfigParameter<double>("min_relative_permeability_liquid");
    auto const min_relative_permeability_gas =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__k_rel_min_gas}
        config.getConfigParameter<double>("min_relative_permeability_gas");
    auto const exponent =
        //! \ogs_file_param{prj__media__medium__properties__property__RelPermBrooksCorey__lambda}
        config.getConfigParameter<double>("lambda");
    if (exponent == 0.)
    {
        OGS_FATAL("Exponent 'lambda' must be positive.");
    }

    return std::make_unique<RelPermBrooksCorey>(
        residual_liquid_saturation,
        residual_gas_saturation,
        min_relative_permeability_liquid,
        min_relative_permeability_gas,
        exponent);
}

}  // namespace MaterialPropertyLib
