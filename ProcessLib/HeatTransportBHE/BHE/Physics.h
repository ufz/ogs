/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE  // namespace of borehole heat exchanger
{
inline double prandtlNumber(double const& viscosity,
                            double const& heat_capacity,
                            double const& heat_conductivity)
{
    return viscosity * heat_capacity / heat_conductivity;
}

inline double reynoldsNumber(double const velocity_norm,
                             double const pipe_diameter,
                             double const viscosity,
                             double const density)
{
    return velocity_norm * pipe_diameter / (viscosity / density);
}

inline double nusseltNumber(double const reynolds_number,
                            double const prandtl_number,
                            double const pipe_diameter,
                            double const pipe_length)
{
    if (reynolds_number < 2300.0)
    {
        return 4.364;
    }
    if (reynolds_number < 10000.0)
    {
        double const gamma = (reynolds_number - 2300) / (10000 - 2300);

        return (1.0 - gamma) * 4.364 +
               gamma *
                   ((0.0308 / 8.0 * 1.0e4 * prandtl_number) /
                    (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) *
                               (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
                    (1.0 + std::pow(pipe_diameter / pipe_length, 2.0 / 3.0)));
    }

    double const xi = std::pow(1.8 * std::log10(reynolds_number) - 1.5, -2.0);
    return (xi / 8.0 * reynolds_number * prandtl_number) /
           (1.0 + 12.7 * std::sqrt(xi / 8.0) *
                      (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
           (1.0 + std::pow(pipe_diameter / pipe_length, 2.0 / 3.0));
}

// Pipe aspect ratio is the pipe's diameter equivalent over pipe's length.
inline double nusseltNumberAnnulus(double const reynolds_number,
                                   double const prandtl_number,
                                   double const diameter_ratio,
                                   double const pipe_aspect_ratio)
{
    if (reynolds_number < 2300.0)
    {
        return 3.66 + (4.0 - 0.102 / (diameter_ratio + 0.02)) *
                          std::pow(diameter_ratio, 0.04);
    }
    if (reynolds_number < 10000.0)
    {
        double const gamma = (reynolds_number - 2300) / (10000 - 2300);

        return (1.0 - gamma) *
                   (3.66 + (4.0 - 0.102 / (diameter_ratio + 0.02))) *
                   std::pow(diameter_ratio, 0.04) +
               gamma *
                   ((0.0308 / 8.0 * 1.0e4 * prandtl_number) /
                    (1.0 + 12.7 * std::sqrt(0.0308 / 8.0) *
                               (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
                    (1.0 + std::pow(pipe_aspect_ratio, 2.0 / 3.0)) *
                    ((0.86 * std::pow(diameter_ratio, 0.84) + 1.0 -
                      0.14 * std::pow(diameter_ratio, 0.6)) /
                     (1.0 + diameter_ratio)));
    }
    double const xi = std::pow(1.8 * std::log10(reynolds_number) - 1.5, -2.0);
    return (xi / 8.0 * reynolds_number * prandtl_number) /
           (1.0 + 12.7 * std::sqrt(xi / 8.0) *
                      (std::pow(prandtl_number, 2.0 / 3.0) - 1.0)) *
           (1.0 + std::pow(pipe_aspect_ratio, 2.0 / 3.0)) *
           ((0.86 * std::pow(diameter_ratio, 0.84) + 1.0 -
             0.14 * std::pow(diameter_ratio, 0.6)) /
            (1.0 + diameter_ratio));
}
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
