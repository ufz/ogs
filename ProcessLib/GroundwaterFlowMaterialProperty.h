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
#include "MaterialLib/Fluid/Viscosity/ViscosityType.h"
#include "MaterialLib/Fluid/Storage/StorageType.h"
#include "MaterialLib/DensityType.h"

namespace ProcessLib
{
using Permeability = MaterialLib::ParameterBase<MaterialLib::PermeabilityType>;
using Storage = MaterialLib::ParameterBase<MaterialLib::StorageType>;
using Density = MaterialLib::ParameterBase<MaterialLib::DensityType>;

struct PorousMediumProperty
{
    /// Hydraulic conductivity for groundwater flow process.
    std::unique_ptr<Permeability> conductivity;

    /// Storage.
    std::unique_ptr<Storage> storage;

    // Porosity (add hereafter if it is needed)
};

/// Class of material parameters of groundwater flow process.
class GroundwaterFlowMaterialProperty
{
        using ConfigTree = boost::property_tree::ptree;
    public:
        GroundwaterFlowMaterialProperty(ConfigTree const& config);

        ~GroundwaterFlowMaterialProperty() {}

        Permeability *getPermeability( const std::size_t mat_id) const
        {
            return _porous_medium_property[mat_id]->conductivity.get();
        }

        Storage *getStorage(std::size_t mat_id) const
        {
            return _porous_medium_property[mat_id]->storage.get();
        }

        Density *getDensity() const
        {
            return _density.get();
        }


    private:
        /// Density.
        std::unique_ptr<Density> _density;

        /// Porous medium properties for each material zone.
        std::vector<std::unique_ptr<PorousMediumProperty>> _porous_medium_property;
};

}   // namespace

#endif  // PROCESS_LIB_GROUNDWATER_FLOW_MATERIAL_PROPERTY_H_
