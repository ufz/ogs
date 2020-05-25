/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#include "VanGenuchtenCapillaryPressureSaturation.h"

#include <algorithm>
#include <cmath>

#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace PorousMedium
{
double VanGenuchtenCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    if (has_regularized_)
    {
        double Sg = 1 - saturation;
        if (Sg <= 1 - saturation_r_ && Sg >= saturation_nonwet_r_)
        {
            return getPcBarvGSg(Sg);
        }
        if (Sg < saturation_nonwet_r_)
        {
            return getPcBarvGSg(saturation_nonwet_r_) +
                   getdPcdSvGBar(saturation_nonwet_r_) *
                       (Sg - saturation_nonwet_r_);
        }

        return getPcBarvGSg(1 - saturation_r_) +
               getdPcdSvGBar(1 - saturation_r_) * (Sg - 1 + saturation_r_);
    }
    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double Se = (S - saturation_r_) / (saturation_max_ - saturation_r_);
    const double pc = pb_ * std::pow(std::pow(Se, (-1.0 / m_)) - 1.0, 1.0 - m_);
    return std::clamp(pc, minor_offset_, pc_max_);
}

double VanGenuchtenCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc = std::max(minor_offset_, capillary_pressure);
    const double Se = std::pow(std::pow(pc / pb_, 1.0 / (1.0 - m_)) + 1.0, -m_);
    const double S = Se * (saturation_max_ - saturation_r_) + saturation_r_;
    return std::clamp(S, saturation_r_ + minor_offset_,
                      saturation_max_ - minor_offset_);
}

double VanGenuchtenCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    if (has_regularized_)
    {
        double const Sg = 1 - saturation;
        if (Sg >= saturation_nonwet_r_ && Sg <= 1 - saturation_r_)
        {
            return -getdPcdSvGBar(Sg);
        }
        if (Sg < saturation_nonwet_r_)
        {
            return -getdPcdSvGBar(saturation_nonwet_r_);
        }

        return -getdPcdSvGBar(1 - saturation_r_);
    }
    if (saturation < saturation_r_)
    {
        return 0;
    }
    if (saturation > saturation_max_)
    {
        return 0;
    }

    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double val1 = std::pow(
        ((S - saturation_r_) / (saturation_max_ - saturation_r_)), -1.0 / m_);
    const double val2 = std::pow(val1 - 1.0, -m_);
    return pb_ * (m_ - 1.0) * val1 * val2 / (m_ * (S - saturation_r_));
}

double VanGenuchtenCapillaryPressureSaturation::getd2PcdS2(
    const double saturation) const
{
    if (has_regularized_)
    {
        OGS_FATAL(
            "Second derivative of regularized van-Genuchten saturation "
            "pressure relation is not implemented.");
    }
    if (saturation < saturation_r_)
    {
        return 0;
    }
    if (saturation > saturation_max_)
    {
        return 0;
    }

    const double S = std::clamp(saturation, saturation_r_ + minor_offset_,
                                saturation_max_ - minor_offset_);
    const double val1 = std::pow(
        ((S - saturation_r_) / (saturation_max_ - saturation_r_)), 1.0 / m_);
    return -pb_ / (m_ * m_ * (S - saturation_r_) * (S - saturation_r_)) *
           std::pow(1 - val1, -m_ - 1) * std::pow(val1, m_ - 1) *
           ((1 - m_ * m_) * val1 + m_ - 1);
}
/// Regularized van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getPcBarvGSg(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::saturation_nonwet_r_;
    double const S_lr = CapillaryPressureSaturation::saturation_r_;
    double const S_bar = getSBar(Sg);
    return getPcvGSg(S_bar) - getPcvGSg(Sg_r + (1 - Sg_r - S_lr) * xi_ / 2);
}
/// Regularized van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getSBar(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::saturation_nonwet_r_;
    double const S_lr = CapillaryPressureSaturation::saturation_r_;
    return Sg_r + (1 - xi_) * (Sg - Sg_r) + 0.5 * xi_ * (1 - Sg_r - S_lr);
}
///  van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getPcvGSg(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::saturation_nonwet_r_;
    double const S_lr = CapillaryPressureSaturation::saturation_r_;
    double const S_le = (1 - Sg - S_lr) /
                        (1 - Sg_r - CapillaryPressureSaturation::saturation_r_);
    return pb_ * std::pow(std::pow(S_le, (-1.0 / m_)) - 1.0, 1.0 - m_);
}
/// derivative dPCdS based on regularized van Genuchten capillary
/// pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getdPcdSvGBar(double Sg) const
{
    double S_bar = getSBar(Sg);
    return getdPcdSvG(S_bar) * (1 - xi_);
}
/// derivative dPCdS based on standard van Genuchten capillary
/// pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getdPcdSvG(
    const double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::saturation_nonwet_r_;
    double const S_lr = CapillaryPressureSaturation::saturation_r_;
    double const nn = 1 / (1 - m_);
    double const S_le = (1 - Sg - S_lr) / (1 - Sg_r - S_lr);
    return pb_ * (1 / (m_ * nn)) * (1 / (1 - S_lr - Sg_r)) *
           std::pow(std::pow(S_le, (-1 / m_)) - 1, (1 / nn) - 1) *
           std::pow(S_le, (-1 / m_)) / S_le;
}

}  // namespace PorousMedium
}  // namespace MaterialLib
