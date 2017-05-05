/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include "RichardsFlowMaterialProperties.h"
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
        std::unique_ptr<RichardsFlowMaterialProperties>&& material_,
        Eigen::VectorXd const  specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        Parameter<double> const& temperature_)
        : material(std::move(material_)),
          specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          temperature(temperature_)
    {
    }

    RichardsFlowProcessData(RichardsFlowProcessData&& other)
        : material(std::move(other.material)),
          specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          temperature(other.temperature)
    {
    }

    //! Copies are forbidden.
    RichardsFlowProcessData(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(RichardsFlowProcessData&&) = delete;

    std::unique_ptr<RichardsFlowMaterialProperties> material;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const has_mass_lumping;
    Parameter<double> const& temperature;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
