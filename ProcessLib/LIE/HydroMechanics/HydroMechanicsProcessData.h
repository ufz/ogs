/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
#define PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_

#include <memory>
#include <Eigen/Dense>

#include "MeshLib/PropertyVector.h"

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "MaterialLib/FractureModels/FractureModelBase.h"

#include "ProcessLib/ProcessVariable.h"
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
template <unsigned GlobalDim>
struct HydroMechanicsProcessData
{
    HydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>&&
            material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& specific_storage_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, GlobalDim, 1> const& specific_body_force_,
        std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>>&& fracture_model,
        std::unique_ptr<FractureProperty>&& fracture_prop,
        Parameter<double> const& initial_effective_stress_,
        Parameter<double> const& initial_fracture_effective_stress_,
        Parameter<double> const& initial_pressure_,
        bool const use_initial_stress_as_reference_
        )
        : material{std::move(material_)},
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_),
          fracture_model{std::move(fracture_model)},
          fracture_property{std::move(fracture_prop)},
          initial_effective_stress(initial_effective_stress_),
          initial_fracture_effective_stress(initial_fracture_effective_stress_),
          initial_pressure(initial_pressure_),
          use_initial_stress_as_reference(use_initial_stress_as_reference_)
    {
    }

    HydroMechanicsProcessData(HydroMechanicsProcessData&& other)
        : material{std::move(other.material)},
          intrinsic_permeability(other.intrinsic_permeability),
          specific_storage(other.specific_storage),
          fluid_viscosity(other.fluid_viscosity),
          fluid_density(other.fluid_density),
          biot_coefficient(other.biot_coefficient),
          porosity(other.porosity),
          solid_density(other.solid_density),
          specific_body_force(other.specific_body_force),
          fracture_model{std::move(other.fracture_model)},
          fracture_property{std::move(other.fracture_property)},
          initial_effective_stress(other.initial_effective_stress),
          initial_fracture_effective_stress(other.initial_fracture_effective_stress),
          initial_pressure(other.initial_pressure),
          use_initial_stress_as_reference(other.use_initial_stress_as_reference),
          dt(other.dt),
          t(other.t)
    {
    }

    //! Copies are forbidden.
    HydroMechanicsProcessData(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HydroMechanicsProcessData&&) = delete;

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<GlobalDim>>
        material;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& specific_storage;
    Parameter<double> const& fluid_viscosity;
    Parameter<double> const& fluid_density;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& porosity;
    Parameter<double> const& solid_density;
    Eigen::Matrix<double, GlobalDim, 1> const specific_body_force;
    std::unique_ptr<MaterialLib::Fracture::FractureModelBase<GlobalDim>> fracture_model;
    std::unique_ptr<FractureProperty> fracture_property;

    Parameter<double> const& initial_effective_stress;
    Parameter<double> const& initial_fracture_effective_stress;
    Parameter<double> const& initial_pressure;
    bool const use_initial_stress_as_reference;

    double dt = 0.0;
    double t = 0.0;
    ProcessLib::ProcessVariable const* pv_p = nullptr;

    // mesh properties for output
    MeshLib::PropertyVector<double>* mesh_prop_stress_xx = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_yy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_stress_xy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xx = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_yy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_strain_xy = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_velocity = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_b = nullptr;
    MeshLib::PropertyVector<double>* mesh_prop_k_f = nullptr;
};

}  // namespace HydroMechanics
}  // namespace LIE
}  // namespace ProcessLib

#endif  // PROCESSLIB_LIE_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
