/*!
  \file  ProcessParameterFluid.cpp
  \brief Define the members of class ProcessParameterFluid in this file.

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

#include "ProcessParameterFluid.h"

namespace ProcessLib
{

ProcessParameterFluid::ProcessParameterFluid(ConfigTree const& config)
    : _name(config.get<std::string>("name"))
{
    DBUG("Set fluid parameters for %s", this->_name.c_str());

    // Permeabiliy
    {
        auto const &perm_config = config.find("permeability");
        auto const &perm_k_config = perm_config  == config.not_found() ?
                                    perm_config : config.find("hydraulic_conductivity");
        if(perm_k_config == config.not_found())
        {
            WARN("No permeability or hydraulic_conductivity found.");
        }
    }
}

}   // namespace
