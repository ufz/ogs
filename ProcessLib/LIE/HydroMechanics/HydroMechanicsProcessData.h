// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>

#include "MaterialLib/FractureModels/FractureModelBase.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MeshLib/ElementStatus.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/LIE/Common/FractureProperty.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace LIE
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids;
    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<DisplacementDim>>
        fracture_model;
    std::vector<FractureProperty> fracture_properties;

    ParameterLib::Parameter<double> const& initial_effective_stress;
    ParameterLib::Parameter<double> const& initial_fracture_effective_stress;

    bool const deactivate_matrix_in_flow;

    /// An indicator to use the B bar method \cite hughes1980generalization to
    /// tackle the  volumetric locking.
    const bool use_b_bar;

    std::vector<JunctionProperty> junction_properties = {};

    MeshLib::PropertyVector<int> const* mesh_prop_materialIDs = nullptr;
    std::vector<int> map_materialID_to_fractureID = {};

    // a table of connected fracture IDs for each element
    std::vector<std::vector<int>> vec_ele_connected_fractureIDs = {};
    std::vector<std::vector<int>> vec_ele_connected_junctionIDs = {};

    std::unique_ptr<MeshLib::ElementStatus> p_element_status = nullptr;
    ParameterLib::Parameter<double> const* p0 = nullptr;

    // mesh properties for output
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* element_velocities = nullptr;
    MeshLib::PropertyVector<double>* element_local_jumps = nullptr;
    MeshLib::PropertyVector<double>* element_fracture_stresses = nullptr;
    MeshLib::PropertyVector<double>* element_fracture_velocities = nullptr;

    MeshLib::PropertyVector<double>* mesh_prop_b = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_k_f = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_fracture_shear_failure = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_p = nullptr;

    MeshLib::PropertyVector<double>* mesh_prop_nodal_forces = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_nodal_forces_jump = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_hydraulic_flow = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
