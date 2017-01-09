/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace MeshLib
{
class Element;
}

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
        std::unique_ptr<MaterialLib::TwoPhaseFlowWithPP::
                            TwoPhaseFlowWithPPMaterialProperties>&& material_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Pc_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_wet_,
        MathLib::PiecewiseLinearInterpolation const& interpolated_Kr_nonwet_)
        : _specific_body_force(specific_body_force_),
          _has_gravity(has_gravity_),
          _has_mass_lumping(has_mass_lumping_),
          _material(std::move(material_)),
          _interpolated_Pc(interpolated_Pc_),
          _interpolated_Kr_wet(interpolated_Kr_wet_),
          _interpolated_Kr_nonwet(interpolated_Kr_nonwet_)
    {
    }

    TwoPhaseFlowWithPPProcessData(TwoPhaseFlowWithPPProcessData&& other)
        : _specific_body_force(other._specific_body_force),
          _has_gravity(other._has_gravity),
          _has_mass_lumping(other._has_mass_lumping),
          _material(std::move(other._material)),
          _interpolated_Pc(other._interpolated_Pc),
          _interpolated_Kr_wet(other._interpolated_Kr_wet),
          _interpolated_Kr_nonwet(other._interpolated_Kr_nonwet)
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
    Eigen::VectorXd const _specific_body_force;

    bool const _has_gravity;

    //! Enables lumping of the mass matrix.
    bool const _has_mass_lumping;
    std::unique_ptr<
        MaterialLib::TwoPhaseFlowWithPP::TwoPhaseFlowWithPPMaterialProperties>
        _material;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Pc;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Kr_wet;
    MathLib::PiecewiseLinearInterpolation const& _interpolated_Kr_nonwet;
};

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
