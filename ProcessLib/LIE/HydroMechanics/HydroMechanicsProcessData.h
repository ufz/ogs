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
template <int GlobalDim>
struct HydroMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids;
    std::map<int,
             std::shared_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>>
        solid_materials;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    Eigen::Matrix<double, GlobalDim, 1> const specific_body_force;
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>
        fracture_model;
    std::unique_ptr<FractureProperty> fracture_property;
    ParameterLib::Parameter<double> const& initial_effective_stress;
    ParameterLib::Parameter<double> const& initial_fracture_effective_stress;

    bool const deactivate_matrix_in_flow;

    /// An indicator to use the B bar method \cite hughes1980generalization to
    /// tackle the  volumetric locking.
    const bool use_b_bar;

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
