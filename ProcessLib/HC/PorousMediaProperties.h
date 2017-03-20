/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include <Eigen/Dense>
#include "BaseLib/reorderVector.h"


#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "ProcessLib/Parameter/SpatialPosition.h"

namespace ProcessLib
{
namespace HC
{

class PorousMediaProperties
{
public:
    PorousMediaProperties(
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
            porosity_models,
        std::vector<Eigen::MatrixXd>&& intrinsic_permeability_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
            specific_storage_models,
        std::vector<int>&& material_ids)
        : _porosity_models(std::move(porosity_models)),
          _intrinsic_permeability_models(
              std::move(intrinsic_permeability_models)),
          _specific_storage_models(std::move(specific_storage_models)),
          _material_ids(std::move(material_ids))
    {
    }

    PorousMediaProperties(PorousMediaProperties&& other)
        : _porosity_models(std::move(other._porosity_models)),
          _intrinsic_permeability_models(other._intrinsic_permeability_models),
          _specific_storage_models(std::move(other._specific_storage_models)),
          _material_ids(other._material_ids)
    {
    }

    MaterialLib::PorousMedium::Porosity const& getPorosity(
        double t, SpatialPosition const& pos) const;

    Eigen::MatrixXd const& getIntrinsicPermeability(
        double t, SpatialPosition const& pos) const;

    MaterialLib::PorousMedium::Storage const& getSpecificStorage(
        double t, SpatialPosition const& pos) const;

private:
    int getMaterialID(SpatialPosition const& pos) const;
private:
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<Eigen::MatrixXd> _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _specific_storage_models;
    std::vector<int> _material_ids;
};

}
}
