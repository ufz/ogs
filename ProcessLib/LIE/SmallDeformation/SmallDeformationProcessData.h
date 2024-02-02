/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "MaterialLib/FractureModels/FractureModelBase.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MeshLib/PropertyVector.h"
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
          _fracture_model{std::move(fracture_model)},
          fracture_properties(std::move(fracture_properties_)),
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
    std::vector<FractureProperty> fracture_properties;
    std::vector<JunctionProperty> junction_properties;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::vector<int> _map_materialID_to_fractureID;

    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> _vec_ele_connected_fractureIDs;
    std::vector<std::vector<int>> _vec_ele_connected_junctionIDs;

    // mesh properties to output element's stress.
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* element_local_jumps = nullptr;
    MeshLib::PropertyVector<double>* element_fracture_stresses = nullptr;

    // mesh property for fracture aperture
    MeshLib::PropertyVector<double>* _mesh_prop_b = nullptr;

    double const _reference_temperature;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
