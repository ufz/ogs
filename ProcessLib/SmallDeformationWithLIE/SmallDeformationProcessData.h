/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONPROCESSDATA_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONPROCESSDATA_H_

#include <memory>

#include "MeshLib/PropertyVector.h"

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/FractureModels/FractureModelBase.h"

#include "ProcessLib/SmallDeformationWithLIE/Common/FractureProperty.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{
template <int DisplacementDim>
struct SmallDeformationProcessData
{
    SmallDeformationProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&& material,
        std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>>&& fracture_model,
        std::unique_ptr<FractureProperty>&& fracture_prop
        )
        : _material{std::move(material)}, _fracture_model{std::move(fracture_model)},
          _fracture_property{std::move(fracture_prop)}
    {
    }

    SmallDeformationProcessData(SmallDeformationProcessData&& other)
        : _material{std::move(other._material)},
          _fracture_model{std::move(other._fracture_model)},
          _fracture_property{std::move(other._fracture_property)}
    {
    }

    //! Copies are forbidden.
    SmallDeformationProcessData(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>> _material;
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>> _fracture_model;
    std::unique_ptr<FractureProperty> _fracture_property;

    double dt = 0.0;
    double t = 0.0;

    MeshLib::PropertyVector<double>* _mesh_prop_stress_xx = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_yy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_xy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_xx = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_yy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_xy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_b = nullptr;
};

}  // namespace SmallDeformationWithLIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_SMALLDEFORMATIONPROCESSDATA_H_
