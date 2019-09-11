/**
 * \author Norbert Grunwald
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
 * \class RelPermLiakopoulos
 * \brief Relative permeability function for Liakopoulos Benchmark
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function
 * of capillary pressure.
 */
class RelPermLiakopoulos final : public Property
{
private:
    Medium* _medium;
    /**
Parameters for Liakopoulos relative permeability:
Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
formulations and coupling strategies for two-phase flow in deforming
porous media. Comput. Geosci., p. 24ff.

Those parameters are fixed for that particular model, no need to change them.
*/
    const double _residual_liquid_saturation = 0.2;
    const double _parameter_a = 2.207;
    const double _parameter_b = 1.0121;
    const double _exponent = 3.;
    const double _min_relative_permeability_gas = 1.0e-4;

public:
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
                "The property 'RelPermLiakopoulos' is implemented on the "
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

}  // namespace MaterialPropertyLib
