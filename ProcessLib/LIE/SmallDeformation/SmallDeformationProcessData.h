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
    MeshLib::PropertyVector<int> const* const material_ids;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>>
        fracture_model;
    std::vector<FractureProperty> fracture_properties;

    double const reference_temperature;

    /// An indicator to use the B bar method \cite hughes1980generalization to
    /// tackle the  volumetric locking.
    const bool use_b_bar;

    std::vector<JunctionProperty> junction_properties = {};

    MeshLib::PropertyVector<int> const* mesh_prop_materialIDs = nullptr;
    std::vector<int> map_materialID_to_fractureID = {};

    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> vec_ele_connected_fractureIDs = {};
    std::vector<std::vector<int>> vec_ele_connected_junctionIDs = {};

    // mesh properties to output element's stress.
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* element_local_jumps = nullptr;
    MeshLib::PropertyVector<double>* element_fracture_stresses = nullptr;

    // mesh property for fracture aperture
    MeshLib::PropertyVector<double>* mesh_prop_b = nullptr;
};

}  // namespace SmallDeformation
}  // namespace LIE
}  // namespace ProcessLib
