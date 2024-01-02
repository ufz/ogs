/**
 * \file
 * \author Norbert Grunwald
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree-fwd.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class RelPermLiakopoulos
 * \brief Relative permeability function for the wetting phase
 *        of the Liakopoulos experiment.
 *
 * \details This property must be a medium property, it
 * computes the permeability reduction due to saturation as function
 * of capillary pressure.
 *
 * The relative permeability is given by an empirical formula
 * \f[k^{\mathrm{rel}}_{\mathrm{L}}=1 - 2.207(1 - s_L)^{1.0121}\f]
 *
 * with effective saturation
 * \f[s_{\mathrm{e}}=\frac{s_{\mathrm{L}}-s^{\mathrm{r}}_{\mathrm{L}}}{1 -
s^{\mathrm{r}}_{\mathrm{L}}}\f]
 *
 * where
 * \f[\lambda=3.0\f]
 * \f[s^{\mathrm{r}}_{\mathrm{L}}=0.2\f]
 */
class RelPermLiakopoulos final : public Property
{
private:
    /**
    Parameters for Liakopoulos relative permeability:
    Asadi, R., Ataie-Ashtiani, B. (2015): A Comparison of finite volume
    formulations and coupling strategies for two-phase flow in deforming
    porous media. Comput. Geosci., p. 24ff.

    Those parameters are fixed for that particular model, no need to change
    them.
    */
    const double residual_liquid_saturation_ = 0.2;
    const double maximal_liquid_saturation_ = 1.;
    const double parameter_a_ = 2.207;
    const double parameter_b_ = 1.0121;

public:
    explicit RelPermLiakopoulos(std::string name);

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'RelPermLiakopoulos' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property values_ and dValues_.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
