/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "TwoPhaseFlowWithPPMaterialProperties.h"

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace TwoPhaseFlowWithPP
{
struct TwoPhaseFlowWithPPProcessData
{
    TwoPhaseFlowWithPPProcessData(
        Eigen::VectorXd const specific_body_force_,
        bool const has_gravity_,
        bool const has_mass_lumping_,
        ParameterLib::Parameter<double> const& temperature_,
        std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>&& material_,
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_)
        : specific_body_force(specific_body_force_),
          has_gravity(has_gravity_),
          has_mass_lumping(has_mass_lumping_),
          temperature(temperature_),
          material(std::move(material_)),
          media_map(std::move(media_map_))

    {
    }

    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData&& other)
        : specific_body_force(other.specific_body_force),
          has_gravity(other.has_gravity),
          has_mass_lumping(other.has_mass_lumping),
          temperature(other.temperature),
          material(std::move(other.material)),
          media_map(std::move(other.media_map))
    {
    }

    //! Copies are forbidden.
    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData const&) =
        delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(TwoPhaseFlowWithPPProcessData&&) = delete;

    //! Specific body forces applied to solid and fluid.
    //! It is usually used to apply gravitational forces.
    //! A vector of displacement dimension's length.
    Eigen::VectorXd const specific_body_force;

    bool const has_gravity;

    //! Enables lumping of the mass matrix.
    bool const has_mass_lumping;
    ParameterLib::Parameter<double> const& temperature;
    std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties> material;
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
