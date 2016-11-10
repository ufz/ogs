/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
#define PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_

#include <Eigen/Dense>

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
struct HydroMechanicsProcessData
{
    HydroMechanicsProcessData(
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>&&
            material_,
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& specific_storage_,
        Parameter<double> const& fluid_viscosity_,
        Parameter<double> const& fluid_density_,
        Parameter<double> const& biot_coefficient_,
        Parameter<double> const& porosity_,
        Parameter<double> const& solid_density_,
        Eigen::Matrix<double, DisplacementDim, 1> const& specific_body_force_)
        : material{std::move(material_)},
          intrinsic_permeability(intrinsic_permeability_),
          specific_storage(specific_storage_),
          fluid_viscosity(fluid_viscosity_),
          fluid_density(fluid_density_),
          biot_coefficient(biot_coefficient_),
          porosity(porosity_),
          solid_density(solid_density_),
          specific_body_force(specific_body_force_)
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

    std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
        material;
    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& specific_storage;
    Parameter<double> const& fluid_viscosity;
    Parameter<double> const& fluid_density;
    Parameter<double> const& biot_coefficient;
    Parameter<double> const& porosity;
    Parameter<double> const& solid_density;
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;
    double dt;
    double t;
};

}  // namespace HydroMechanics
}  // namespace ProcessLib

#endif  // PROCESSLIB_HYDROMECHANICS_HYDROMECHANICSPROCESSDATA_H_
