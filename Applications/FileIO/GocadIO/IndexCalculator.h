/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
        : _x_dim(x_dim),
          _y_dim(y_dim),
          _z_dim(z_dim),
          _n_nodes(x_dim * y_dim * z_dim),
          _n_cells((_x_dim - 1) * (_y_dim - 1) * (_z_dim - 1))
    {
    }

    IndexCalculator() = default;

    std::size_t operator()(std::array<std::size_t, 3> const& c) const
    {
        const std::size_t idx(c[2] * _x_dim * _y_dim + c[1] * _x_dim + c[0]);
        if (idx >= _n_nodes)
        {
            return std::numeric_limits<std::size_t>::max();
        }
        return idx;
    }

    std::size_t getCellIdx(std::size_t u, std::size_t v, std::size_t w) const
    {
        // ensure (u,v,w) is a valid cell
        if (u >= _x_dim - 1 || v >= _y_dim - 1 || w >= _z_dim - 1)
        {
            ERR("GocadSGridReader::IndexCalculator::getCellIdx(): At least "
                "one grid coordinate to big.");
            ERR("\t Given: ({:d}, {:d}, {:d}), max allowed cell grid coords: "
                "({:d}, {:d}, {:d}).",
                u, v, w, _x_dim - 1, _y_dim - 1, _z_dim - 1);
            return std::numeric_limits<std::size_t>::max();
        }

        return (_x_dim - 1) * (_y_dim - 1) + v * (_x_dim - 1) + u;
    }

    std::array<std::size_t, 3> getCoordsForID(std::size_t id) const
    {
        std::array<std::size_t, 3> const coords{
            (id % (_x_dim * _y_dim)) % _x_dim,
            (id % (_x_dim * _y_dim)) / _x_dim, id / (_x_dim * _y_dim)};
        return coords;
    }

    std::size_t _x_dim{0};
    std::size_t _y_dim{0};
    std::size_t _z_dim{0};
    std::size_t _n_nodes{0};
    std::size_t _n_cells{0};
};

}  // end namespace Gocad
}  // end namespace FileIO
