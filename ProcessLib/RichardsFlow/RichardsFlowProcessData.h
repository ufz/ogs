/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_RICHARDSFLOW_RICHARDSFLOWPROCESSDATA_H
#define PROCESSLIB_RICHARDSFLOW_RICHARDSFLOWPROCESSDATA_H

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace RichardsFlow
{
struct RichardsFlowProcessData
{
    RichardsFlowProcessData(
        Parameter<double> const& intrinsic_permeability_,
        Parameter<double> const& porosity_,
        Parameter<double> const& viscosity_,
        Parameter<double> const& storage_,
        Parameter<double> const& water_density_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        std::map<std::string,
                 std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
            curves_)
        : intrinsic_permeability(intrinsic_permeability_),
          porosity(porosity_),
          viscosity(viscosity_),
          storage(storage_),
          water_density(water_density_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          curves(curves_)
    {
    }

    RichardsFlowProcessData(RichardsFlowProcessData&& other)
        : intrinsic_permeability(other.intrinsic_permeability),
          porosity(other.porosity),
          viscosity(other.viscosity),
          storage(other.storage),
          water_density(other.water_density),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          curves(other.curves)
    {
    }

    //! Copies are forbidden.
    RichardsFlowProcessData(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData&&) = delete;

    Parameter<double> const& intrinsic_permeability;
    Parameter<double> const& porosity;
    Parameter<double> const& viscosity;
    Parameter<double> const& storage;
    Parameter<double> const& water_density;
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves;
    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib

#endif  // PROCESSLIB_RichardsFlow_RichardsFlowPROCESSDATA_H
