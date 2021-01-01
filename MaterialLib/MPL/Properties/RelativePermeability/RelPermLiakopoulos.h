/**
 * \file
 * \author Norbert Grunwald
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
 *
 * Wetting (liquid) phase relative permeability is given by the empirical
 * formula
 * \f[k^{\mathrm{rel}}_{\mathrm{L}}=1 - 2.207(1 - s_L)^{1.0121}\f]
 *
 * Non-wetting (gas) phase permeability is computed using the Brooks-Corey
 * model
 * \f[k^{\mathrm{rel}}_{\mathrm{G}}=(1 -
s_{\mathrm{e}})^2\left(1-s_{\mathrm{e}}^{\frac{\left(2+\lambda\right)}{\lambda}}\right)\f]
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
    const double exponent_ = 3.;
    const double min_relative_permeability_gas_ = 1.0e-4;
    const double dse_dsL_ =
        1. / (maximal_liquid_saturation_ - residual_liquid_saturation_);

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
