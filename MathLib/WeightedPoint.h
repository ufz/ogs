/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <array>
#include <iosfwd>
#include <limits>

namespace MathLib
{
//! Represents a point of a certain dimension that has a weight attached.
//!
//! Used, e.g., in numerical quadrature.
class WeightedPoint
{
public:
    //! Constructs a 0...3D weighted point depending on the passed coordinates
    //! array.
    template <std::size_t dim>
    WeightedPoint(std::array<double, dim> const& coords, double const weight)
        : weight_{weight}, dim_{dim}
    {
        static_assert(dim <= 3);
        std::size_t i = 0;
        if constexpr (dim > 0)  // avoids compiler warning
        {
            for (; i < dim; ++i)
            {
                coords_[i] = coords[i];
            }
        }
        for (; i < 3; ++i)
        {
            // fill the rest with NaN for safety reasons
            coords_[i] = std::numeric_limits<double>::quiet_NaN();
        }
    }

    //! Constructs a 0D weighted point.
    explicit WeightedPoint(double const weight) : weight_{weight}, dim_{0}
    {
        // fill with NaN for safety reasons
        coords_.fill(std::numeric_limits<double>::quiet_NaN());
    }

    bool operator==(WeightedPoint const& other) const
    {
        if (weight_ != other.weight_)
        {
            return false;
        }
        if (dim_ != other.dim_)
        {
            return false;
        }
        for (std::size_t comp = 0; comp < dim_; ++comp)
        {
            if (coords_[comp] != other.coords_[comp])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!=(WeightedPoint const& other) const
    {
        return !(*this == other);
    }

    const double* data() const { return coords_.data(); }

    double getWeight() const { return weight_; }

    //! The point dimension, i.e., the number of its coordinates.
    std::size_t getDimension() const { return dim_; }

    //! Access a specific coordinate.
    double operator[](std::size_t coord_idx) const
    {
        return coords_[coord_idx];
    }

private:
    double weight_;
    std::array<double, 3> coords_;
    std::size_t dim_;
};

std::ostream& operator<<(std::ostream& os, MathLib::WeightedPoint const& wp);

}  // namespace MathLib
