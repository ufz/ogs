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
    HydroMechanicsProcessData(
        MeshLib::PropertyVector<int> const* const material_ids_,
        std::map<
            int,
            std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>>&&
            solid_materials_,
        ParameterLib::Parameter<double> const& intrinsic_permeability_,
        ParameterLib::Parameter<double> const& specific_storage_,
        ParameterLib::Parameter<double> const& fluid_viscosity_,
        ParameterLib::Parameter<double> const& fluid_density_,
        ParameterLib::Parameter<double> const& biot_coefficient_,
        ParameterLib::Parameter<double> const& porosity_,
        ParameterLib::Parameter<double> const& solid_density_,
        Eigen::Matrix<double, GlobalDim, 1>
            specific_body_force_,
        std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>&&
            fracture_model,
        std::unique_ptr<FracturePropertyHM>&& fracture_prop,
        ParameterLib::Parameter<double> const& initial_effective_stress_,
        ParameterLib::Parameter<double> const&
            initial_fracture_effective_stress_,
        bool const deactivate_matrix_in_flow_,
        double const reference_temperature_)
        : material_ids(material_ids_),
          solid_materials{std::move(solid_materials_)},
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          specific_body_force(std::move(specific_body_force_)),
          fracture_model{std::move(fracture_model)},
          fracture_property{std::move(fracture_prop)},
          initial_effective_stress(initial_effective_stress_),
          initial_fracture_effective_stress(initial_fracture_effective_stress_),
          deactivate_matrix_in_flow(deactivate_matrix_in_flow_),
          reference_temperature(reference_temperature_)
    {
    }

    HydroMechanicsProcessData(HydroMechanicsProcessData&& other) = default;

    //! Copies are forbidden.
    HydroMechanicsProcessData(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData&&) = delete;

    MeshLib::PropertyVector<int> const* const material_ids;
    std::map<int,
             std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& intrinsic_permeability;
    ParameterLib::Parameter<double> const& specific_storage;
    ParameterLib::Parameter<double> const& fluid_viscosity;
    ParameterLib::Parameter<double> const& fluid_density;
    ParameterLib::Parameter<double> const& biot_coefficient;
    ParameterLib::Parameter<double> const& porosity;
    ParameterLib::Parameter<double> const& solid_density;
    Eigen::Matrix<double, GlobalDim, 1> const specific_body_force;
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>
        fracture_model;
    std::unique_ptr<FracturePropertyHM> fracture_property;
    ParameterLib::Parameter<double> const& initial_effective_stress;
    ParameterLib::Parameter<double> const& initial_fracture_effective_stress;

    bool const deactivate_matrix_in_flow;
    std::unique_ptr<MeshLib::ElementStatus> p_element_status;
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

    double const reference_temperature;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib
