/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>

#include "MaterialLib/PorousMedium/Permeability/Permeability.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace PorousMedium
{
/// The purpose of the class DupuitPermeability is the implementation of the
/// special diffusion coefficient \f$h K(h)\f$, where \f$h\f$ is the hydraulic
/// head of the previous non-linear iteration and \f$K\f$ is the hydraulic
/// conductivity. The diffusion coefficient is used in unconfined groundwater
/// flow equation.
class DupuitPermeability final : public Permeability
{
public:
    DupuitPermeability(
        ParameterLib::Parameter<double> const& permeability_parameter,
        int const dimension)
        : Permeability(permeability_parameter, dimension)
    {
    }

    /// @copydoc Permeability::getValue()
    ///
    /// The third parameter variable is (mis)used as aquifer thickness.
    Eigen::MatrixXd getValue(const double t,
                             ParameterLib::SpatialPosition const& pos,
                             const double variable,
                             const double temperature) const override
    {
        return variable * Permeability::getValue(t, pos, variable, temperature);
    }
};

}  // namespace PorousMedium
}  // namespace MaterialLib
