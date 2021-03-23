/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    : length_(L), dL0_(dL0), max_dL_(max_dL), multiplier_(multiplier)
{
    // Check if accumulated subdivisions can ever sum up to length.
    // Cf. geometric series formula.
    if (multiplier < 1.0 && dL0 / (1.0 - multiplier) < L) {
        OGS_FATAL(
            "Using dL0={:g} and multiplier={:g} the generated subdivisions can "
            "not sum up to a total length of {:g}.",
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
        x += std::min(max_dL_,
                      dL0_ * std::pow(multiplier_, static_cast<double>(i)));
        i++;
    } while (x < length_);

    if (vec_x.back() < length_)
    {
        double last_dx = vec_x[vec_x.size() - 1] - vec_x[vec_x.size() - 2];
        if (length_ - vec_x.back() < last_dx)
        {
            vec_x[vec_x.size() - 1] = length_;
        }
        else
        {
            vec_x.push_back(length_);
        }
    }
    return vec_x;
}

GradualSubdivisionFixedNum::GradualSubdivisionFixedNum(
    const double L, const std::size_t num_subdivisions, const double multiplier)
    : length_{L}, num_subdivisions_{num_subdivisions}, multiplier_{multiplier}
{
}

std::vector<double> GradualSubdivisionFixedNum::operator()() const
{
    std::vector<double> subdivisions;
    subdivisions.reserve(num_subdivisions_ + 1);
    subdivisions.push_back(0.0);
    auto const q = multiplier_;

    if (q == 1.0) {
        double const dx = length_ / num_subdivisions_;

        for (std::size_t i = 1; i < num_subdivisions_; ++i)
        {
            subdivisions.push_back(dx * i);
        }
    } else {
        // compute initial subdivision size
        auto const a =
            length_ * (q - 1.0) / (std::pow(q, num_subdivisions_) - 1.0);

        double qi = q;  // q^i
        for (std::size_t i = 1; i < num_subdivisions_; ++i)
        {
            subdivisions.push_back(a * (qi - 1.0) / (q - 1.0));
            qi *= q;
        }
    }

    subdivisions.push_back(length_);

    return subdivisions;
}

}  // namespace BaseLib
