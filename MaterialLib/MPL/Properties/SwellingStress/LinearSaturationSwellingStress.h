/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 8:19 AM
 */

#pragma once

#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
class Phase;
/**
 * \brief This class defines a linear saturation dependent swelling stress
 *   model for the materials that swell strongly when water content increases.
 *
 *   Clay materials like bentonite have a high swelling capacity in dry state,
 * and their swelling property can be described by this model.
 *
 *   The model takes the form
 *    \f[  {\mathbf{\sigma}}^{\text{sw}} =
 * {\alpha}_{\text{sw}}  (S-S_0) \mathbf{I} \f]
 *  where
 *  \f${\alpha}_{\text{sw}}\f$ is a coefficient,  and \f$S_0\f$ is the
 *  initial saturation.
 *
 *   The coefficient gives the swelling stress at full saturation, which can be
 *  computed as
 *    \f[ {\alpha}_{\text{sw}} =
 * \frac{{{\sigma}}^{\text{sw}}_{\text{max}}}{(S_{\text{max}}-S_0)}
 *     \f]
 *  where \f${{\sigma}}^{\text{sw}}_{\text{max}}\f$  represents the swelling
 * stress at full saturation, and \f$S_{\text{max}}\f$ is the maximum
 * saturation.
 *
 *
 *  In the numerical analysis, the stress always takes the incremental form.
 *  Therefore the model becomes
 *   \f[\Delta  {\mathbf{\sigma}}^{\text{sw}} =
 * {\alpha}_{\text{sw}}  \Delta S \mathbf{I} \f]
 * with only one parameter.
 *
 */
class LinearSaturationSwellingStress final : public Property
{
public:
    LinearSaturationSwellingStress(std::string name, double const coefficient);

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'LinearSaturationSwellingStress' is "
                "implemented on the 'phase' scale only.");
        }
    }

    /// \return The swelling stress
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// \return The derivative of the swelling stress with respect to
    /// saturation.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /** The coefficient of the model, which gives the swelling stress at full
     saturation.
     *  By this way, the coefficient can be obtained by
     * \f[ {\alpha}_{\text{sw}} =
     * \frac{{{\sigma}}^{\text{sw}}_{\text{max}}}{(S_{\text{max}}-S_0)}
       \f]
     * with  \f${{\sigma}}^{\text{sw}}_{\text{max}}\f$   the swelling
     * stress at full saturation, and \f$S_{\text{max}}\f$ is the maximum
     saturation.
     */
    double const coefficient_;
};

}  // namespace MaterialPropertyLib
