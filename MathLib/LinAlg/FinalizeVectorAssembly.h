/**
 * \file FinalizeVectorAssembly.h
 * \author Wenqing Wang
 * \date Oct, 2013
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FINALIZEVECTORASSEMBLY_H_
#define FINALIZEVECTORASSEMBLY_H_

namespace MathLib
{
/// General function to finalize the vector assembly
template <typename VEC_T>
void finalizeVectorAssembly(VEC_T &)
{
}

} // MathLib

#endif
