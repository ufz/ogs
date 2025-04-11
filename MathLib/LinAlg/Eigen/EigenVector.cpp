/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenVector.h"

#include <range/v3/algorithm/copy.hpp>

#include "BaseLib/Error.h"
#include "EigenMapTools.h"

#ifndef NDEBUG
#include <fstream>
#endif

namespace MathLib
{
void EigenVector::copyValues(std::vector<double>& u) const
{
    u.resize(size());
    toVector(u) = vec_;
}

void EigenVector::copyValues(std::span<double> u) const
{
    if (u.size() != static_cast<std::size_t>(vec_.size()))
    {
        OGS_FATAL(
            "EigenVector copy values; size mismatch. Trying to copy a vector "
            "of size {:d} to a span of size {:d}.",
            size(), u.size());
    }
    ranges::copy(vec_, u.begin());
}

#ifndef NDEBUG
void EigenVector::write(const std::string& filename) const
{
    std::ofstream os(filename);
    os << vec_;
}
#endif

}  // namespace MathLib
