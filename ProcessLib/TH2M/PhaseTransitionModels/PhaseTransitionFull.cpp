/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionFull.h"

#include "MaterialLib/PhysicalConstant.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransitionFull::PhaseTransitionFull(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModel(media)
{
    DBUG("Create PhaseTransitionFull constitutive model.");
}

PhaseTransitionModelVariables PhaseTransitionFull::updateConstitutiveVariables(
    PhaseTransitionModelVariables const& /*phase_transition_model_variables*/,
    const MaterialPropertyLib::Medium* /*medium*/,
    MaterialPropertyLib::VariableArray /*variables*/,
    ParameterLib::SpatialPosition /*pos*/, double const /*t*/,
    const double /*dt*/) const
{
    OGS_FATAL(
        "PhaseTransitionFull::updateConstitutiveVariables is not implemented.");
}

}  // namespace TH2M
}  // namespace ProcessLib
