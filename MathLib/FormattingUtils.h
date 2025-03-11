/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <spdlog/fmt/bundled/ostream.h>
#include <spdlog/fmt/bundled/ranges.h>

#include <Eigen/Core>
#include <concepts>

#include "mathlib_export.h"

namespace MathLib
{
struct EigenIOFormat
{
    static MATHLIB_EXPORT const Eigen::IOFormat full_precision;
};
}  // namespace MathLib

// disable fmt's range formatting for Eigen types
template <typename T>
    requires std::derived_from<T, Eigen::DenseBase<T>>
struct fmt::is_range<T, char> : std::false_type
{
};

template <typename T>
    requires std::derived_from<T, Eigen::DenseBase<T>>
struct fmt::formatter<T> : fmt::ostream_formatter
{
    auto format(T const& value, fmt::format_context& ctx) const
    {
        return fmt::ostream_formatter::format(
            value.format(MathLib::EigenIOFormat::full_precision), ctx);
    }
};
