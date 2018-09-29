/**
 * @file FinalizeMatrixAssembly.h
 * @author Thomas Fischer
 * @date Jun 11, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

namespace MathLib
{

template <typename MAT_T>
bool finalizeMatrixAssembly(MAT_T &)
{
    return true;
}

} // MathLib
