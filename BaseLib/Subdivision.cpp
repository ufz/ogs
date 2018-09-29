/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Subdivision.h"

#include <algorithm>
#include <cmath>

#include <BaseLib/Error.h>

namespace BaseLib
{
GradualSubdivision::GradualSubdivision(const double L,
                                       const double dL0,
                                       const double max_dL,
                                       const double multiplier)
    : _length(L), _dL0(dL0), _max_dL(max_dL), _multiplier(multiplier)
{
    // Check if accumulated subdivisions can ever sum up to length.
    // Cf. geometric series formula.
    if (multiplier < 1.0 && dL0 / (1.0 - multiplier) < L) {
        OGS_FATAL(
            "Using dL0=%g and multiplier=%g the generated subdivisions can not "
            "sum up to a total length of %g.",
            dL0,
            multiplier,
            L);
    }
}

std::vector<double> GradualSubdivision::operator()() const
{
    std::vector<double> vec_x;

    double x = 0;
    unsigned i = 0;
    do {
        vec_x.push_back(x);
        x += std::min(_max_dL,
                      _dL0 * std::pow(_multiplier, static_cast<double>(i)));
        i++;
    } while (x < _length);

    if (vec_x.back() < _length) {
        double last_dx = vec_x[vec_x.size() - 1] - vec_x[vec_x.size() - 2];
        if (_length - vec_x.back() < last_dx)
            vec_x[vec_x.size() - 1] = _length;
        else
            vec_x.push_back(_length);
    }
    return vec_x;
}

GradualSubdivisionFixedNum::GradualSubdivisionFixedNum(
    const double L, const std::size_t num_subdivisions, const double multiplier)
    : _length{L}, _num_subdivisions{num_subdivisions}, _multiplier{multiplier}
{
}

std::vector<double> GradualSubdivisionFixedNum::operator()() const
{
    std::vector<double> subdivisions;
    subdivisions.reserve(_num_subdivisions + 1);
    subdivisions.push_back(0.0);
    auto const q = _multiplier;

    if (q == 1.0) {
        double const dx = _length / _num_subdivisions;

        for (std::size_t i = 1; i < _num_subdivisions; ++i) {
            subdivisions.push_back(dx * i);
        }
    } else {
        // compute initial subdivision size
        auto const a =
            _length * (q - 1.0) / (std::pow(q, _num_subdivisions) - 1.0);

        double qi = q;  // q^i
        for (std::size_t i = 1; i < _num_subdivisions; ++i) {
            subdivisions.push_back(a * (qi - 1.0) / (q - 1.0));
            qi *= q;
        }
    }

    subdivisions.push_back(_length);

    return subdivisions;
}

}  // namespace BaseLib
