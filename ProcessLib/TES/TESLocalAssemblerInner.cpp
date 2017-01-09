/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements
 * from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#include "TESLocalAssemblerInner-fwd.h"

#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES
#include "TESLocalAssemblerInner-impl.h"
#endif

namespace ProcessLib
{
namespace TES
{
#ifdef OGS_EIGEN_DYNAMIC_SHAPE_MATRICES

template class TESLocalAssemblerInner<
    LocalAssemblerTraits<ShapeMatrixPolicyType<void, 0>, 0, 0, 0>>;

static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 1,
              "inconsistent use of macros");
#else
static_assert(OGS_EIGEN_DYNAMIC_SHAPE_MATRICES_FLAG == 0,
              "inconsistent use of macros");
#endif
}
}
