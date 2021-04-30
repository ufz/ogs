/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PhaseTransitionDissolution.h"

namespace ProcessLib
{
namespace TH2M
{
PhaseTransitionDissolution::PhaseTransitionDissolution(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media)
    : PhaseTransitionModels(media)
{
    DBUG("Create PhaseTransitionDissolution constitutive model.");
}

void PhaseTransitionDissolution::computeConstitutiveVariables(
    const MaterialPropertyLib::Medium* /*medium*/,
    MaterialPropertyLib::VariableArray /*variables*/,
    ParameterLib::SpatialPosition /*pos*/, double const /*t*/,
    const double /*dt*/)
{
    OGS_FATAL(
        "PhaseTransitionDissolution::computeConstitutiveVariables is not "
        "implemented.");
}

}  // namespace TH2M
}  // namespace ProcessLib
