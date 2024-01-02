/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Core>

#include "NumericalStabilization.h"

namespace NumLib
{
namespace detail
{

inline Eigen::MatrixXd getHydrodynamicDispersion(
    Eigen::MatrixXd const& pore_diffusion_coefficient,
    Eigen::VectorXd const& velocity,
    double const porosity,
    double const solute_dispersivity_transverse,
    double const solute_dispersivity_longitudinal)
{
    double const velocity_magnitude = velocity.norm();
    if (velocity_magnitude == 0.0)
    {
        return porosity * pore_diffusion_coefficient;
    }

    auto const dim = velocity.size();
    Eigen::MatrixXd const& I(Eigen::MatrixXd::Identity(dim, dim));
    return porosity * pore_diffusion_coefficient +
           solute_dispersivity_transverse * velocity_magnitude * I +
           (solute_dispersivity_longitudinal - solute_dispersivity_transverse) /
               velocity_magnitude * velocity * velocity.transpose();
}

inline Eigen::MatrixXd getHydrodynamicDispersionWithArtificialDiffusion(
    IsotropicDiffusionStabilization const& stabilizer,
    std::size_t const element_id,
    Eigen::MatrixXd const& pore_diffusion_coefficient,
    Eigen::VectorXd const& velocity,
    double const porosity,
    double const solute_dispersivity_transverse,
    double const solute_dispersivity_longitudinal)
{
    double const velocity_magnitude = velocity.norm();
    if (velocity_magnitude == 0.0)
    {
        return porosity * pore_diffusion_coefficient;
    }

    double const artificial_diffusion =
        stabilizer.computeArtificialDiffusion(element_id, velocity_magnitude);

    auto const dim = velocity.size();
    Eigen::MatrixXd const& I(Eigen::MatrixXd::Identity(dim, dim));
    return porosity * pore_diffusion_coefficient +
           (solute_dispersivity_transverse * velocity_magnitude +
            artificial_diffusion) *
               I +
           (solute_dispersivity_longitudinal - solute_dispersivity_transverse) /
               velocity_magnitude * velocity * velocity.transpose();
}
}  // namespace detail

inline Eigen::MatrixXd computeHydrodynamicDispersion(
    NumericalStabilization const& stabilizer,
    std::size_t const element_id,
    Eigen::MatrixXd const& pore_diffusion_coefficient,
    Eigen::VectorXd const& velocity,
    double const porosity,
    double const solute_dispersivity_transverse,
    double const solute_dispersivity_longitudinal)
{
    return std::visit(
        [&](auto&& stabilizer)
        {
            using Stabilizer = std::decay_t<decltype(stabilizer)>;
            if constexpr (std::is_same_v<Stabilizer,
                                         IsotropicDiffusionStabilization>)
            {
                return detail::getHydrodynamicDispersionWithArtificialDiffusion(
                    stabilizer,
                    element_id,
                    pore_diffusion_coefficient,
                    velocity,
                    porosity,
                    solute_dispersivity_transverse,
                    solute_dispersivity_longitudinal);
            }

            return detail::getHydrodynamicDispersion(
                pore_diffusion_coefficient,
                velocity,
                porosity,
                solute_dispersivity_transverse,
                solute_dispersivity_longitudinal);
        },
        stabilizer);
}

}  // namespace NumLib
