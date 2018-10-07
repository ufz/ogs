/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct PipeParameters
{
    /**
     * radius of the pipline inner side
     * unit is m
     */
    double const r_inner;

    /**
     * radius of the pipline outer side
     * unit is m
     */
    double const r_outer;

    /**
     * pipe-in wall thickness
     * unit is m
     */
    double const b_in;

    /**
     * pipe-out wall thickness
     * unit is m
     */
    double const b_out;

    /**
     * thermal conductivity of the pipe wall
     * unit is kg m sec^-3 K^-1
     */

    double const lambda_p;

    /**
     * thermal conductivity of the inner pipe wall
     * unit is kg m sec^-3 K^-1
     */
    double const lambda_p_i;

    /**
     * thermal conductivity of the outer pipe wall
     * unit is kg m sec^-3 K^-1
     */
    double const lambda_p_o;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
