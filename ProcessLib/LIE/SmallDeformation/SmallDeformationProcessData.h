/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MeshLib/PropertyVector.h"

#include "MaterialLib/FractureModels/FractureModelBase.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"

#include "ProcessLib/LIE/Common/FractureProperty.h"
#include "ProcessLib/LIE/Common/JunctionProperty.h"

namespace MeshLib
{
class Element;
}  // namespace MeshLib

namespace ProcessLib
{
namespace LIE
{
namespace SmallDeformation
{
template <int DisplacementDim>
struct SmallDeformationProcessData
{
    SmallDeformationProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<int,
                 std::unique_ptr<
                     MaterialLib::Solids::MechanicsBase<DisplacementDim>>>&&
            solid_materials_,
        std::unique_ptr<
            MaterialLib::Fracture::FractureModelBase<DisplacementDim>>&&
            fracture_model,
        std::vector<FractureProperty>&& fracture_properties_,
        double const reference_temperature)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          fracture_model_{std::move(fracture_model)},
          fracture_properties(std::move(fracture_properties_)),
          reference_temperature_(reference_temperature)
    {
    }

    SmallDeformationProcessData(SmallDeformationProcessData&& other) = default;

    //! Copies are forbidden.
    SmallDeformationProcessData(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;

    /// The constitutive relation for the mechanical part.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>>
        fracture_model_;
    std::vector<FractureProperty> fracture_properties;
    std::vector<JunctionProperty> junction_properties;

    MeshLib::PropertyVector<int> const* mesh_prop_materialIDs_ = nullptr;
    std::vector<int> map_materialID_to_fractureID_;

    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> vec_ele_connected_fractureIDs_;
    std::vector<std::vector<int>> vec_ele_connected_junctionIDs_;

    // mesh properties to output element's stress.
    MeshLib::PropertyVector<double>* mesh_prop_stress_xx_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_yy_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_zz_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_xy_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_xz_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_yz_ = nullptr;

    // mesh properties to output element's strain
    MeshLib::PropertyVector<double>* mesh_prop_strain_xx_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_yy_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_zz_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xy_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xz_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_yz_ = nullptr;

    // mesh property for fracture aperture
    MeshLib::PropertyVector<double>* mesh_prop_b_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_w_n_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_w_s_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_stress_shear_ = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_stress_normal_ =
        nullptr;

    double const reference_temperature_;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
