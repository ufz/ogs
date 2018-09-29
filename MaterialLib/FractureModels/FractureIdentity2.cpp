/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FractureIdentity2.h"

namespace MaterialLib
{
namespace Fracture
{

namespace
{
template <typename VectorType, int DisplacementDim>
VectorType initIdentity2()
{
    VectorType ivec = VectorType::Zero(DisplacementDim);
    ivec[DisplacementDim-1] = 1;
    return ivec;
}
} // anonymous namespace

template <int DisplacementDim>
const typename FractureIdentity2<DisplacementDim>::VectorType
    FractureIdentity2<DisplacementDim>::value = initIdentity2<
        typename FractureIdentity2<DisplacementDim>::VectorType, DisplacementDim>();

// template instantiation
template struct FractureIdentity2<2>;
template struct FractureIdentity2<3>;

}  // namespace Fracture
}  // namespace MaterialLib

