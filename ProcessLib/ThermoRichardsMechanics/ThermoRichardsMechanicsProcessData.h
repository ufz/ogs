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

#include "BaseLib/StrongType.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Common/HydroMechanics/InitialStress.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
/// A type helping to avoid confusion of different boolean values.
using InitializePorosityFromMediumProperty =
    BaseLib::StrongType<bool, struct InitializePorosityFromMediumPropertyTag>;

template <int DisplacementDim, typename ConstitutiveTraits>
struct ThermoRichardsMechanicsProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// The constitutive relation for the mechanical part.
    std::map<int, std::unique_ptr<
                      typename ConstitutiveTraits::SolidConstitutiveRelation>>
        solid_materials;

    InitialStress const initial_stress;

    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    bool const apply_mass_lumping;

    const bool use_TaylorHood_elements;

    bool const apply_body_force_for_deformation;

    InitializePorosityFromMediumProperty const
        initialize_porosity_from_medium_property;

    MeshLib::PropertyVector<double>* element_saturation = nullptr;
    MeshLib::PropertyVector<double>* element_porosity = nullptr;
    MeshLib::PropertyVector<double>* element_liquid_density = nullptr;
    MeshLib::PropertyVector<double>* element_viscosity = nullptr;
    MeshLib::PropertyVector<double>* element_stresses = nullptr;
    MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;
    MeshLib::PropertyVector<double>* pressure_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib
