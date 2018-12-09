/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
        std::vector<FractureProperty>&& vec_fracture_prop,
        double const reference_temperature)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          _fracture_model{std::move(fracture_model)},
          _vec_fracture_property(std::move(vec_fracture_prop)),
          _reference_temperature(reference_temperature)
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
        _fracture_model;
    std::vector<FractureProperty> _vec_fracture_property;
    std::vector<JunctionProperty> _vec_junction_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::vector<int> _map_materialID_to_fractureID;

    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> _vec_ele_connected_fractureIDs;
    std::vector<std::vector<int>> _vec_ele_connected_junctionIDs;

    double dt = 0.0;
    double t = 0.0;

    // mesh properties to output element's stress.
    MeshLib::PropertyVector<double>* _mesh_prop_stress_xx = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_yy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_zz = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_xy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_xz = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_stress_yz = nullptr;

    // mesh properties to output element's strain
    MeshLib::PropertyVector<double>* _mesh_prop_strain_xx = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_yy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_zz = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_xy = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_xz = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_strain_yz = nullptr;

    // mesh property for fracture aperture
    MeshLib::PropertyVector<double>* _mesh_prop_b = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_w_n = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_w_s = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_fracture_stress_shear = nullptr;
    MeshLib::PropertyVector<double>* _mesh_prop_fracture_stress_normal =
        nullptr;

    double const _reference_temperature;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
