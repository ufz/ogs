/*!
   \file  GroundwaterFlowMaterialProperty.h
   \brief Declare a class of material parameters of groundwater flow process.

   \author Wenqing Wang
   \date Feb 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef PROCESS_LIB_GROUNDWATER_FLOW_MATERIAL_PROPERTY_H_
#define PROCESS_LIB_GROUNDWATER_FLOW_MATERIAL_PROPERTY_H_

#include <boost/property_tree/ptree.hpp>

#include <vector>
#include <memory>

#include "MaterialLib/ParameterBase.h"
#include "MaterialLib/Fluid/Permeability/PermeabilityType.h"
#include "MaterialLib/Fluid/Storage/StorageType.h"
#include "MaterialLib/DensityType.h"

namespace ProcessLib
{
using Permeability = MaterialLib::ParameterBase<MaterialLib::PermeabilityType>;
using Storage = MaterialLib::ParameterBase<MaterialLib::StorageType>;
using Density = MaterialLib::ParameterBase<MaterialLib::DensityType>;

/// Class of material parameters of groundwater flow process.
class GroundwaterFlowMaterialProperty
{
        using ConfigTree = boost::property_tree::ptree;
    public:
        GroundwaterFlowMaterialProperty(ConfigTree const& config);

        ~GroundwaterFlowMaterialProperty();

        Permeability *getConductivity( const std::size_t mat_id) const
        {
            return _conductivity[mat_id];
        }

        Storage *getStorage(std::size_t mat_id) const
        {
            return _storage[mat_id];
        }

        Density *getDensity() const
        {
            return _density.get();
        }


    private:
        /// Density.
        std::unique_ptr<Density> _density;

        /// Hydraulic conductivity for groundwater flow process.
        std::vector<Permeability*> _conductivity;

        /// Storage.
        std::vector<Storage*> _storage;

        // Porosity (add hereafter if it is needed)
};

}   // namespace

#endif  // PROCESS_LIB_GROUNDWATER_FLOW_MATERIAL_PROPERTY_H_
