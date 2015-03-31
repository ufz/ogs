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

#include <algorithm> // for std::sort

#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp> // for boost::split
#include "logog/include/logog.hpp"

#include "GroundwaterFlowMaterialProperty.h"

namespace ProcessLib
{

GroundwaterFlowMaterialProperty::GroundwaterFlowMaterialProperty(ConfigTree const& config)
{
    // Porous medium properties
    auto const& pmp_config = config.find("porous_medium_property");
    if (pmp_config == config.not_found())
        WARN("No porous medium property found.");

    _conductivity.resize(pmp_config->second.size());

    for (auto const& pmp_iterator : pmp_config->second)
    {
        ConfigTree const& pmp_config = pmp_iterator.second;
        std::size_t const mat_id = pmp_config.get<std::size_t>("id");

        auto const &hydraulic_conductivity_config = pmp_config.get_child("hydraulic_conductivity");
        const std::string hy_type = hydraulic_conductivity_config.get<std::string>("type");
        if(hy_type.find("isotropic") != std::string::npos && hy_type.find("aniso") == std::string::npos)
        {
            const double k = hydraulic_conductivity_config.get<double>("value");

            using IsotropicHydraulicConductivity =
                MaterialLib::TensorParameter<MaterialLib::PermeabilityType,
                MaterialLib::ConstantScalarModel, double>;
            IsotropicHydraulicConductivity *K
                = new IsotropicHydraulicConductivity(k);

            _conductivity[mat_id] = K;
        }
        else if(hy_type.find("anisotropic") != std::string::npos)
        {
            const std::string values_in_str = hydraulic_conductivity_config.get<std::string>("values");
            std::vector<std::string> strs_tok;
            boost::split(strs_tok, values_in_str, boost::is_any_of(" "));
            std::vector<double> values;
            for(auto str : strs_tok) // blank spaces between data are allowed.
            {
                if(str.empty())
                    continue;
                values.push_back(std::stod(str));
            }

            const std::size_t data_size= values.size();
            if(data_size != 4 || data_size != 9)
            {
                ERR("Number of values for an anisotropic hydraulic conductivity tensor must be 4 or 9.");
            }
            const int dim = (data_size == 4) ? 2 : 3;
            Matrix K(dim, dim);
            for(int i=0; i<dim; i++)
            {
                for(int j=0; j<dim; j++)
                {
                    K(i, j) = values[i*dim + j];
                }
            }
            using AnisotropicHydraulicConductivity
            = MaterialLib::TensorParameter<MaterialLib::PermeabilityType,
            MaterialLib::ConstantTensor<Matrix>, Matrix>;
            AnisotropicHydraulicConductivity *anisK = new AnisotropicHydraulicConductivity(K);

            _conductivity[mat_id] = anisK;
        }
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
