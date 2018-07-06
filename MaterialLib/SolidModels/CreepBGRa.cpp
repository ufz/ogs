/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   CreepBGRa.cpp
 *  Created on July 6, 2018, 9:53 AM
 */

#include "CreepBGRa.h"

#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace Solids
{
namespace Creep
{

template <int DisplacementDim>
boost::optional<std::tuple<typename CreepBGRa<DisplacementDim>::KelvinVector,
                           std::unique_ptr<typename MechanicsBase<
                               DisplacementDim>::MaterialStateVariables>,
                           typename CreepBGRa<DisplacementDim>::KelvinMatrix>>
CreepBGRa<DisplacementDim>::integrateStress(
    double const t,
    ProcessLib::SpatialPosition const& x,
    double const dt,
    KelvinVector const& eps_prev,
    KelvinVector const& eps,
    KelvinVector const& sigma_prev,
    typename MechanicsBase<DisplacementDim>::MaterialStateVariables const&
        material_state_variables, double const T, double const p)
{
    return LinearElasticIsotropic<DisplacementDim>::integrateStress(
        t, x, dt, eps_prev, eps, sigma_prev, material_state_variables, T, p);
}

template class CreepBGRa<2>;
template class CreepBGRa<3>;

}  // namespace Creep
}  // namespace Solids
}  // namespace MaterialLib
