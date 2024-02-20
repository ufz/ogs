/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <spdlog/fmt/bundled/ostream.h>

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
