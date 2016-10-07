/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   LiquidFlowMaterialProperties.cpp
 *
 * Created on August 18, 2016, 11:49 AM
 */

#include "LiquidFlowMaterialProperties.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Mesh.h"

#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

namespace ProcessLib
{
namespace LiquidFlow
{
LiquidFlowMaterialProperties::LiquidFlowMaterialProperties(
    BaseLib::ConfigTree const& config)
{
    DBUG("Reading material properties of liquid flow process.");

    //! \ogs_file_param{prj__material_property__fluid}
    auto const& fluid_config = config.getConfigSubtree("fluid");

    // Get fluid properties
    //! \ogs_file_param{prj__material_property__fluid__density}
    auto const& rho_conf = fluid_config.getConfigSubtree("density");
    density_l = MaterialLib::Fluid::createFluidDensityModel(rho_conf);
    //! \ogs_file_param{prj__material_property__fluid__viscosity}
    auto const& mu_conf = fluid_config.getConfigSubtree("viscosity");
    viscosity = MaterialLib::Fluid::createViscosityModel(mu_conf);

    // Get porous properties
    //! \ogs_file_param{prj__material_property__porous_medium}
    auto const& poro_config = config.getConfigSubtree("porous_medium");
    //! \ogs_file_param{prj__material_property__porous_medium__porous_medium}
    for (auto const& conf : poro_config.getConfigSubtreeList("porous_medium"))
    {
        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__permeability}
        auto const& perm_conf = conf.getConfigSubtree("permeability");
        intrinsic_permeabiliy.emplace_back(
            MaterialLib::PorousMedium::createPermeabilityModel(perm_conf));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__porosity}
        auto const& poro_conf = conf.getConfigSubtree("porosity");
        auto n = MaterialLib::PorousMedium::createPorosityModel(poro_conf);
        porosity.emplace_back(std::move(n));

        //! \ogs_file_param{prj__material_property__porous_medium__porous_medium__storage}
        auto const& stora_conf = conf.getConfigSubtree("storage");
        auto beta = MaterialLib::PorousMedium::createStorageModel(stora_conf);
        storage.emplace_back(std::move(beta));
    }
}

double LiquidFlowMaterialProperties::getMassCoefficient(
    const double porosity_variable, const double storage_variable,
    const double p, const double T, const unsigned material_group_id) const
{
    ArrayType vars;
    vars[0] = T;
    vars[1] = p;
    return porosity[material_group_id]->getValue(porosity_variable, T) *
               density_l->getdValue(
                   vars, MaterialLib::Fluid::PropertyVariableType::pl)
               // Divided by rho_l because the PDE is scaled with rho_l
               / density_l->getValue(vars) +
           storage[material_group_id]->getValue(storage_variable);
}

}  // end of namespace
}  // end of namespace
