/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenVector.h"

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

#ifndef NDEBUG
void EigenVector::write(const std::string& filename) const
{
    std::ofstream os(filename);
    os << vec_;
}
#endif

}  // namespace MathLib
