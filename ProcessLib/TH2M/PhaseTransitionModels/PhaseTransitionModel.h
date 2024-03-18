/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>

#include "../ConstitutiveRelations/PhaseTransitionData.h"
#include "MaterialLib/MPL/Medium.h"

namespace ProcessLib
{
namespace TH2M
{
struct PhaseTransitionModel
{
    explicit PhaseTransitionModel(
        std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const&
            media)
    {
        DBUG("Create phase transition models...");

        // check for minimum requirement definitions in media object
        std::array const required_gas_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};
        std::array const required_liquid_properties = {
            MaterialPropertyLib::viscosity, MaterialPropertyLib::density};

        for (auto const& m : media)
        {
            checkRequiredProperties(m.second->phase("Gas"),
                                    required_gas_properties);
            checkRequiredProperties(m.second->phase("AqueousLiquid"),
                                    required_liquid_properties);
        }
    }

    virtual ~PhaseTransitionModel() = default;

    virtual void updateConstitutiveVariables(
        ConstitutiveRelations::PhaseTransitionData& cv,
        const MaterialPropertyLib::Medium* medium,
        MaterialPropertyLib::VariableArray variables,
        ParameterLib::SpatialPosition pos, double const t,
        double const dt) const = 0;
};

int numberOfComponents(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name);

int findComponentIndex(
    std::map<int, std::shared_ptr<MaterialPropertyLib::Medium>> const& media,
    std::string phase_name, MaterialPropertyLib::PropertyType property_type);

}  // namespace TH2M
}  // namespace ProcessLib
