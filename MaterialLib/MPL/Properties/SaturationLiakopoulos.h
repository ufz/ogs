/**
 * \author Norbert Grunwald
 * \date   27.06.2018
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

#include <limits>
#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class SaturationLiakopoulos
 * \brief A well known soil characteristics function
 * \details This property must be a medium property, it
 * computes the saturation of the wetting phase as function
 * of capillary pressure.
 */
class SaturationLiakopoulos final : public Property
{
private:
    Medium* _medium;
    /**
  Parameters for Liakopoulos saturation curve taken from:
  Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
  formulations and coupling strategies for two-phase flow in deforming
  porous media. Comput. Geosci., p. 24ff.

  Those parameters are fixed for that particular model, no need to change them.
*/
    const double _residual_liquid_saturation = 0.2;
    const double _parameter_a = 1.9722e-11;
    const double _parameter_b = 2.4279;
    const double _p_cap_max = std::pow(
        (1. - _residual_liquid_saturation) / _parameter_a, (1. / _parameter_b));

public:
    /// This method assigns a pointer to the material object that is the owner
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
                "The property 'SaturationLiakopoulos' is implemented on the "
                "'media' scale only.");
        }
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

}  // namespace MaterialPropertyLib
