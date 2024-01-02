/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
/**
 * Hydraulic aperture equals the mechanical aperture s.t. multiplication of the
 * permeability by the mechanical aperture yields the cubic law.
 * The volumetric flux of incompressible fluid flow in a laminar regime
 * within a fracture consisting of two parallel and smooth surfaces is given
 * by the cubic law:
 *
 * \f[ Q = W/(u L) b^3/12 \Delta p,\f]
 * where, \f$Q\f$ is the volumetric flux normal to the flow direction, \f$W\f$
 * is the fracture width, \f$u\f$ is the fluid viscosity, \f$L\f$ is the
 * fracture length, \f$p\f$ is the fluid pressure, and \f$b\f$ is the fracture
 * hydraulic aperture.
 *
 * In the context of Darcy’s law, the fracture permeability
 * becomes: \f[ k_f = b^2/12,\f].
 *
 * Ref: He, Xupeng, et al. "A corrected cubic law
 * for single-phase laminar flow through rough-walled fractures." Advances in
 * Water resources 154 (2021): 103984.
 */
struct CubicLawPermeability final : public Property
{
    explicit CubicLawPermeability(std::string name,
                                  ParameterLib::Parameter<double> const& b)
        : _b(b)
    {
        name_ = std::move(name);
    }

    void checkScale() const override
    {
        if (!std::holds_alternative<Medium*>(scale_))
        {
            OGS_FATAL(
                "The property 'CubicLawPermeability' is implemented on the "
                "'media' scaleonly.");
        }
    }

    /// \return The permeability value.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    /// \return The derivative of permeability with respect to the aperture.
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;

private:
    /// fracture aperture
    ParameterLib::Parameter<double> const& _b;
};
}  // namespace MaterialPropertyLib
