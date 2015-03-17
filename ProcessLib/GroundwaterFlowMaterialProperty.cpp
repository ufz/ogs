/*!
  \file  GroundwaterFlowMaterialProperty.cpp
  \brief Define the members of class GroundwaterFlowMaterialProperty in this file.

  \author Wenqing Wang
  \date Feb 2015

  \copyright
   Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
              Distributed under a Modified BSD License.
              See accompanying file LICENSE.txt or
              http://www.opengeosys.org/project/license
*/

#include <boost/property_tree/ptree.hpp>
#include "logog/include/logog.hpp"

#include "GroundwaterFlowMaterialProperty.h"

#include "MaterialLib/ConstantScalarModel.h"
#include "MaterialLib/TensorParameter.h"
#include "MaterialLib/Fluid/Permeability/PermeabilityType.h"

namespace ProcessLib
{
using namespace MaterialLib;

GroundwaterFlowMaterialProperty::GroundwaterFlowMaterialProperty(ConfigTree const& config)
{
    // Porous medium properties
    auto const& pmp_config = config.find("porous_medium_property");
    if (pmp_config == config.not_found())
        WARN("No porous medium property found.");

    std::vector<std::size_t> mat_ids;
    for (auto const& pmp_iterator : pmp_config->second)
    {
        ConfigTree const& pmp_config = pmp_iterator.second;
        std::size_t const mat_id = pmp_config.get<std::size_t>("id");
        mat_ids.push_back(mat_id);

        auto const &hydraulic_conductivity_config = pmp_config.get_child("hydraulic_conductivity");
        std::string const hy_type = hydraulic_conductivity_config.get<std::string>("type");
        if(hy_type.find("isotropic") != std::string::npos)
        {
            const double k = hydraulic_conductivity_config.get<double>("value");

            MaterialLib::TensorParameter<PermeabilityType, MaterialLib::ConstantScalarModel, double> *K
			   = new MaterialLib::TensorParameter<PermeabilityType, MaterialLib::ConstantScalarModel, double>(k);

            _conductivity.push_back(K);
        }

        // sort _conductivity
    }
}

GroundwaterFlowMaterialProperty::~GroundwaterFlowMaterialProperty()
{
    for(auto cd : _conductivity)
        delete cd;

    for(auto st: _storage)
        delete st;
}
}   // namespace
