/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   VanGenuchtenCapillaryPressureSaturation.cpp
 *
 *  Created on October 28, 2016, 6:05 PM
 */

#include "VanGenuchtenCapillaryPressureSaturation.h"

#include <cmath>

#include "MathLib/MathTools.h"
#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace PorousMedium
{
double VanGenuchtenCapillaryPressureSaturation::getCapillaryPressure(
    const double saturation) const
{
    if (_has_regularized)
    {
        double Sg = 1 - saturation;
        if (Sg <= 1 - _saturation_r && Sg >= _saturation_nonwet_r)
        {
            return getPcBarvGSg(Sg);
        }
        if (Sg < _saturation_nonwet_r)
        {
            return getPcBarvGSg(_saturation_nonwet_r) +
                   getdPcdSvGBar(_saturation_nonwet_r) *
                       (Sg - _saturation_nonwet_r);
        }

        return getPcBarvGSg(1 - _saturation_r) +
               getdPcdSvGBar(1 - _saturation_r) * (Sg - 1 + _saturation_r);
    }
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double Se = (S - _saturation_r) / (_saturation_max - _saturation_r);
    const double pc = _pb * std::pow(std::pow(Se, (-1.0 / _m)) - 1.0, 1.0 - _m);
    return MathLib::limitValueInInterval(pc, _minor_offset, _pc_max);
}

double VanGenuchtenCapillaryPressureSaturation::getSaturation(
    const double capillary_pressure) const
{
    const double pc =
        (capillary_pressure < 0.0) ? _minor_offset : capillary_pressure;
    double Se = std::pow(pc / _pb, 1.0 / (1.0 - _m)) + 1.0;
    Se = std::pow(Se, -_m);
    const double S = Se * (_saturation_max - _saturation_r) + _saturation_r;
    return MathLib::limitValueInInterval(S, _saturation_r + _minor_offset,
                                         _saturation_max - _minor_offset);
}

double VanGenuchtenCapillaryPressureSaturation::getdPcdS(
    const double saturation) const
{
    if (_has_regularized)
    {
        double const Sg = 1 - saturation;
        if (Sg >= _saturation_nonwet_r && Sg <= 1 - _saturation_r)
        {
            return -getdPcdSvGBar(Sg);
        }
        if (Sg < _saturation_nonwet_r)
        {
            return -getdPcdSvGBar(_saturation_nonwet_r);
        }

        return -getdPcdSvGBar(1 - _saturation_r);
    }
    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double val1 = std::pow(
        ((S - _saturation_r) / (_saturation_max - _saturation_r)), -1.0 / _m);
    const double val2 = std::pow(val1 - 1.0, -_m);
    return _pb * (_m - 1.0) * val1 * val2 / (_m * (S - _saturation_r));
}

double VanGenuchtenCapillaryPressureSaturation::getd2PcdS2(
    const double saturation) const
{
    if (_has_regularized)
    {
        OGS_FATAL(
            "Second derivative of regularized van-Genuchten saturation "
            "pressure relation is not implemented.");
    }
    if (saturation < _saturation_r)
        return 0;
    if (saturation > _saturation_max)
        return 0;

    const double S =
        MathLib::limitValueInInterval(saturation, _saturation_r + _minor_offset,
                                      _saturation_max - _minor_offset);
    const double val1 = std::pow(
        ((S - _saturation_r) / (_saturation_max - _saturation_r)), 1.0 / _m);
    return -_pb / (_m * _m * (S - _saturation_r) * (S - _saturation_r)) *
           std::pow(1 - val1, -_m - 1) * std::pow(val1, _m - 1) *
           ((1 - _m * _m) * val1 + _m - 1);
}
/// Regularized van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getPcBarvGSg(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::_saturation_nonwet_r;
    double const S_lr = CapillaryPressureSaturation::_saturation_r;
    double const S_bar = getSBar(Sg);
    return getPcvGSg(S_bar) - getPcvGSg(Sg_r + (1 - Sg_r - S_lr) * _xi / 2);
}
/// Regularized van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getSBar(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::_saturation_nonwet_r;
    double const S_lr = CapillaryPressureSaturation::_saturation_r;
    return Sg_r + (1 - _xi) * (Sg - Sg_r) + 0.5 * _xi * (1 - Sg_r - S_lr);
}
///  van Genuchten capillary pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getPcvGSg(double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::_saturation_nonwet_r;
    double const S_lr = CapillaryPressureSaturation::_saturation_r;
    double const S_le = (1 - Sg - S_lr) /
                        (1 - Sg_r - CapillaryPressureSaturation::_saturation_r);
    return _pb * std::pow(std::pow(S_le, (-1.0 / _m)) - 1.0, 1.0 - _m);
}
/// derivative dPCdS based on regularized van Genuchten capillary
/// pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getdPcdSvGBar(double Sg) const
{
    double S_bar = getSBar(Sg);
    return getdPcdSvG(S_bar) * (1 - _xi);
}
/// derivative dPCdS based on standard van Genuchten capillary
/// pressure-saturation Model
double VanGenuchtenCapillaryPressureSaturation::getdPcdSvG(
    const double Sg) const
{
    double const Sg_r = CapillaryPressureSaturation::_saturation_nonwet_r;
    double const S_lr = CapillaryPressureSaturation::_saturation_r;
    double const nn = 1 / (1 - _m);
    double const S_le = (1 - Sg - S_lr) / (1 - Sg_r - S_lr);
    return _pb * (1 / (_m * nn)) * (1 / (1 - S_lr - Sg_r)) *
           std::pow(std::pow(S_le, (-1 / _m)) - 1, (1 / nn) - 1) *
           std::pow(S_le, (-1 / _m)) / S_le;
}

}  // end namespace
}  // end namespace
