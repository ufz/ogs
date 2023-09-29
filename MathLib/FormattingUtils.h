/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <spdlog/fmt/bundled/ostream.h>

#include <Eigen/Core>
#include <concepts>

template <typename T>
requires std::derived_from<T, Eigen::DenseBase<T>>
struct fmt::formatter<T> : fmt::ostream_formatter
{
};
