/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <spdlog/spdlog.h>

#include <array>
#include <cstddef>
#include <limits>
#include <numeric>

#include "BaseLib/Logging.h"

namespace FileIO
{
namespace Gocad
{
/**
 * Class for calculating the index to given 3d position within the
 * structured grid.
 */
class IndexCalculator final
{
public:
    /**
     * Constructor initializes the dimensions.
     * @param x_dim
     * @param y_dim
     * @param z_dim
     */
    IndexCalculator(std::size_t x_dim, std::size_t y_dim, std::size_t z_dim)
        : x_dim_(x_dim),
          y_dim_(y_dim),
          z_dim_(z_dim),
          n_nodes_(x_dim * y_dim * z_dim),
          n_cells_((x_dim_ - 1) * (y_dim_ - 1) * (z_dim_ - 1))
    {
    }

    IndexCalculator() = default;

    std::size_t operator()(std::array<std::size_t, 3> const& c) const
    {
        const std::size_t idx(c[2] * x_dim_ * y_dim_ + c[1] * x_dim_ + c[0]);
        if (idx >= n_nodes_)
        {
            return std::numeric_limits<std::size_t>::max();
        }
        return idx;
    }

    std::size_t getCellIdx(std::size_t u, std::size_t v, std::size_t w) const
    {
        // ensure (u,v,w) is a valid cell
        if (u >= x_dim_ - 1 || v >= y_dim_ - 1 || w >= z_dim_ - 1)
        {
            ERR("GocadSGridReader::IndexCalculator::getCellIdx(): At least "
                "one grid coordinate to big.");
            ERR("\t Given: ({:d}, {:d}, {:d}), max allowed cell grid coords: "
                "({:d}, {:d}, {:d}).",
                u, v, w, x_dim_ - 1, y_dim_ - 1, z_dim_ - 1);
            return std::numeric_limits<std::size_t>::max();
        }

        return (x_dim_ - 1) * (y_dim_ - 1) + v * (x_dim_ - 1) + u;
    }

    std::array<std::size_t, 3> getCoordsForID(std::size_t id) const
    {
        std::array<std::size_t, 3> const coords{
            (id % (x_dim_ * y_dim_)) % x_dim_,
            (id % (x_dim_ * y_dim_)) / x_dim_, id / (x_dim_ * y_dim_)};
        return coords;
    }

    std::size_t x_dim_{0};
    std::size_t y_dim_{0};
    std::size_t z_dim_{0};
    std::size_t n_nodes_{0};
    std::size_t n_cells_{0};
};

}  // end namespace Gocad
}  // end namespace FileIO
