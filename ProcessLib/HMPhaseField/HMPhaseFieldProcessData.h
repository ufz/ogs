/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <memory>
#include <utility>
#ifdef USE_PETSC
#include <petscsys.h>
#endif

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/SolidModels/PhaseFieldBase.h"
#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace HMPhaseField
{
inline void showEnergyAndWork(const double t, double& _elastic_energy,
                              double& _surface_energy, double& _pressure_work)
{
#ifdef USE_PETSC
    double const elastic_energy = _elastic_energy;
    MPI_Allreduce(&elastic_energy, &_elastic_energy, 1, MPI_DOUBLE, MPI_SUM,
                  PETSC_COMM_WORLD);
    double const surface_energy = _surface_energy;
    MPI_Allreduce(&surface_energy, &_surface_energy, 1, MPI_DOUBLE, MPI_SUM,
                  PETSC_COMM_WORLD);
    double const pressure_work = _pressure_work;
    MPI_Allreduce(&pressure_work, &_pressure_work, 1, MPI_DOUBLE, MPI_SUM,
                  PETSC_COMM_WORLD);
#endif

    INFO("Elastic energy: {} Surface energy: {} Pressure work: {} at time: {} ",
         _elastic_energy, _surface_energy, _pressure_work, t);
};
template <int DisplacementDim>
struct HMPhaseFieldProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    std::map<int, std::shared_ptr<
                      MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& residual_stiffness;
    ParameterLib::Parameter<double> const& crack_resistance;
    ParameterLib::Parameter<double> const& crack_length_scale;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_fracture_direction;
    double irreversible_threshold;
    MaterialLib::Solids::Phasefield::PhaseFieldModel phasefield_model;
    MaterialLib::Solids::Phasefield::EnergySplitModel energy_split_model;
    MaterialLib::Solids::Phasefield::SofteningCurve softening_curve;
    double characteristic_length;
    std::unique_ptr<MaterialLib::Solids::Phasefield::DegradationDerivative>
        degradation_derivative;
    double const diffused_range_parameter;
    double const fluid_compressibility;
    double const fracture_threshold;
    double const fracture_permeability_parameter;
    double const fixed_stress_stabilization_parameter;
    double const spatial_stabilization_parameter;
    ParameterLib::Parameter<double> const& width_init;

    MeshLib::PropertyVector<double>* ele_d = nullptr;
    MeshLib::PropertyVector<double>* width = nullptr;

    double elastic_energy = 0.0;
    double surface_energy = 0.0;
    double pressure_work = 0.0;
    int const _phasefield_process_id = 0;
    int const _hydro_process_id = 1;
    int const _mechanics_related_process_id = 2;
};

}  // namespace HMPhaseField
}  // namespace ProcessLib
