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
 * \brief This class defines a linear saturation rate dependent swelling stress
 *   model for the materials that swell strongly when water content increases.
 *
 *   Clay materials like bentonite have a high swelling capacity in dry state,
 * and their swelling property can be described by this model.
 *
 *  The original model was proposed in \cite rutqvist2011implementation
 *  (equations (39) and (40) on pages 758--759). With a simplification of the
 *  parameters of the original formula and introducing a constraint to avoid
 *  shrinkage stress when saturation is below the initial saturation, the model
 *  takes the form
 *    \f[  {\mathbf{\sigma}}^{\text{sw}} =
 * {\alpha}_{\text{sw}}  (S-S_0) \mathbf{I}, \, \forall S \in
 * [S_0, S_\text{max}] \f]
 *  where
 *  \f${\alpha}_{\text{sw}}\f$ is a coefficient, and \f$S_0\f$ is the
 *  initial saturation, and \f$S_{\text{max}}\f$ is the maximum saturation.
 *  The coefficient gives the swelling stress at full saturation, which can be
 *  computed as
 *    \f[ {\alpha}_{\text{sw}} =
 * \frac{{{\sigma}}^{\text{sw}}_{\text{max}}}{(S_{\text{max}}-S_0)}
 *     \f]
 *  where \f${{\sigma}}^{\text{sw}}_{\text{max}}\f$  represents the swelling
 * stress at full saturation.
 *
 *  In the numerical analysis, the stress always takes the incremental form.
 *  Therefore the model becomes as
 *   \f[\Delta  {\mathbf{\sigma}}^{\text{sw}} =
 * {\alpha}_{\text{sw}}  \Delta S \mathbf{I}, \, \forall S \in
 * [S_0, S_\text{max}] \f]
 *
 *  <b>Note</b>:
 *  <ul>
 *   <li>In the property, saturation means water saturation.</li>
 *   <li>The upper limit of saturation is not guaranteed, but
 *       it is required to be less or equal to \f$S_\text{max}\f$
 *       when calling this property. Therefore it is not checked again in this
 *       class.Thus, this model only needs two input parameters:
 *       \f${\alpha}_{\text{sw}}\f$ and \f$S_0\f$.
 *    </li>
 *   </ul>
 *
 **/
class LinearSaturationSwellingStress final : public Property
{
public:
    LinearSaturationSwellingStress(std::string name,
                                   double const coefficient,
                                   double const reference_saturation);

    void checkScale() const override
    {
        if (!std::holds_alternative<Phase*>(scale_))
        {
            OGS_FATAL(
                "The property 'LinearSaturationSwellingStress' is "
                "implemented on the 'phase' scale only.");
        }
    }

    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// \return \f$\Delta  {{\sigma}}^{\text{sw}} \f$.
    PropertyDataType value(VariableArray const& variable_array,
                           VariableArray const& variable_array_prev,
                           ParameterLib::SpatialPosition const& pos,
                           double const t,
                           double const dt) const override;

    /// \return \f$ \frac{ \partial \Delta {{\sigma}}^{\text{sw}}}{\partial
    /// S}\f$.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /** The coefficient of the model, which gives the swelling stress at full
     saturation.
     *  The coefficient can be obtained by
     * \f[ {\alpha}_{\text{sw}} =
     * \frac{{{\sigma}}^{\text{sw}}_{\text{max}}}{(S_{\text{max}}-S_0)}
       \f]
     * with  \f${{\sigma}}^{\text{sw}}_{\text{max}}\f$   the swelling
     * stress at full saturation, and \f$S_{\text{max}}\f$ the maximum
     saturation.
     */
    double const coefficient_;

    /// The reference saturation, at which the swelling stress is
    /// zero.
    double const reference_saturation_;
};

}  // namespace MaterialPropertyLib
