/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TESLocalAssemblerInner.h"

namespace ProcessLib
{
namespace TES
{
#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES

extern template class TESLocalAssemblerInner<
    LocalAssemblerTraits<ShapeMatrixPolicyType<void, 0>, 0, 0, 0>>;

static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 1,
              "inconsistent use of macros");
#else
static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 0,
              "inconsistent use of macros");
#endif
}
}
