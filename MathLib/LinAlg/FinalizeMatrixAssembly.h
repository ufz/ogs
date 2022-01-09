/**
 * \file
 * \author Thomas Fischer
 * \date Jun 11, 2013
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace MathLib
{
template <typename MAT_T>
bool finalizeMatrixAssembly(MAT_T& /*unused*/)
{
    return true;
}

} // MathLib
