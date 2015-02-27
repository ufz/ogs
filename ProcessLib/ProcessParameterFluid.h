/*!
   \file  ProcessParameterFluid.h
   \brief Declare a class of fluid parameters in this file.

   \author Wenqing Wang
   \date Feb 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PROCESS_LIB_PROCESS_PARAMETER_FLUID_H_
#define PROCESS_LIB_PROCESS_PARAMETER_FLUID_H_

#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <memory>

#include "MaterialLib/ParameterBase.h"
#include "MaterialLib/Fluid/Permeability/PermeabilityType.h"
#include "MaterialLib/Fluid/Viscosity/ViscosityType.h"
#include "MaterialLib/Fluid/Storage/StorageType.h"

namespace ProcessLib
{
using Permeability = MaterialLib::ParameterBase<MaterialLib::PermeabilityType>;
using Viscosity = MaterialLib::ParameterBase<MaterialLib::ViscosityType>;
using Storage = MaterialLib::ParameterBase<MaterialLib::StorageType>;

/// Class to keep fluid parameters of all material groups in a simulated domain.
class ProcessParameterFluid
{
        using ConfigTree = boost::property_tree::ptree;
    public:
        ProcessParameterFluid(ConfigTree const& config);

        ~ProcessParameterFluid();

        /// Get process name.
        std::string const& getName() const;

        const std::vector<std::unique_ptr<Permeability>> &getPermeability() const
        {
            return _permeability;
        }

        const std::vector<std::unique_ptr<Viscosity>> &getViscosity() const
        {
            return _viscosity;
        }

        const std::vector<std::unique_ptr<Storage>> &getStorage() const
        {
            return _storage;
        }

    private:
        /// Process name.
        std::string const _name;

        /// Intrisinc permeability, equivalent to conductivity for groundwater flow process.
        std::vector<std::unique_ptr<Permeability>> _permeability;

        /// Viscosity.
        std::vector<std::unique_ptr<Viscosity>> _viscosity;

        /// Storage.
        std::vector<std::unique_ptr<Storage>> _storage;
};

}   // namespace

#endif  // PROCESS_LIB_PROCESS_PARAMETER_FLUID_H_
