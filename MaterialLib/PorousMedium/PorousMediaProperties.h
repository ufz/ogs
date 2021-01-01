/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>

#include "Permeability/Permeability.h"
#include "Porosity/Porosity.h"
#include "Storage/Storage.h"

#include "ParameterLib/SpatialPosition.h"

namespace MaterialLib
{
namespace PorousMedium
{

class PorousMediaProperties
{
public:
    PorousMediaProperties(
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
            porosity_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>&&
            intrinsic_permeability_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
            specific_storage_models,
        MeshLib::PropertyVector<int> const* const material_ids)
        : _porosity_models(std::move(porosity_models)),
          _intrinsic_permeability_models(
              std::move(intrinsic_permeability_models)),
          _specific_storage_models(std::move(specific_storage_models)),
          _material_ids(material_ids)
    {
    }

    PorousMediaProperties(PorousMediaProperties&& other) = default;

    MaterialLib::PorousMedium::Porosity const& getPorosity(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::Permeability const& getIntrinsicPermeability(
        double t, ParameterLib::SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::Storage const& getSpecificStorage(
        double t, ParameterLib::SpatialPosition const& pos) const;

private:
    int getMaterialID(ParameterLib::SpatialPosition const& pos) const;

private:
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _specific_storage_models;
    MeshLib::PropertyVector<int> const* const _material_ids;
};

}  // namespace PorousMedium
}  // namespace MaterialLib
