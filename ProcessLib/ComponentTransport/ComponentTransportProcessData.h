/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "ChemicalProcessData.h"
#include "MaterialLib/MPL/CreateMaterialSpatialDistributionMap.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"

namespace MaterialPropertyLib
{
class Medium;
}

namespace ProcessLib
{
namespace ComponentTransport
{
struct ComponentTransportProcessData
{
    ComponentTransportProcessData(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        Eigen::VectorXd const& specific_body_force_, bool const has_gravity_,
        bool const non_advective_form_,
        std::unique_ptr<ChemicalProcessData>&& chemical_process_data_)
        : media_map(std::move(media_map_)),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          non_advective_form(non_advective_form_),
          chemical_process_data(std::move(chemical_process_data_))
    {
    }

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const non_advective_form;
    std::unique_ptr<ChemicalProcessData> chemical_process_data;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
