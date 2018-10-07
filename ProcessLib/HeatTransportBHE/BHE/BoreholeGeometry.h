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
struct BoreholeGeometry
{
    /**
     * length/depth of the BHE
     * unit is m
     */
    double const length;

    /**
     * diameter of the BHE
     * unit is m
     */
    double const diameter;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
