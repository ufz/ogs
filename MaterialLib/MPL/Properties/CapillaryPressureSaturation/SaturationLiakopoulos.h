/**
 * \file
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
 *
 * Wetting (liquid) phase saturation is given by the empirical relation
 * \f[s^\mathrm{r}_\mathrm{L}=1 - 1.9722\cdot10^{-11}p_\mathrm{cap}^{2.4279}\f]
 *
 */
class SaturationLiakopoulos final : public Property
{
private:
    /**
      Parameters for Liakopoulos saturation curve taken from:
      Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
      formulations and coupling strategies for two-phase flow in deforming
      porous media. Comput. Geosci., p. 24ff.

      Those parameters are fixed for that particular model, no need to change
      them.
    */
    const double residual_liquid_saturation_ = 0.2;
    const double parameter_a_ = 1.9722e-11;
    const double parameter_b_ = 2.4279;
    const double p_cap_max_ = std::pow(
        (1. - residual_liquid_saturation_) / parameter_a_, (1. / parameter_b_));

public:
    explicit SaturationLiakopoulos(std::string name)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'SaturationLiakopoulos' is implemented on the "
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
};

}  // namespace MaterialPropertyLib
