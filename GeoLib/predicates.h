/**
 * \file
 * \author Lars Bilke
 * \date   2025-02-06
 * \brief  Definitions of the predicates functions.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

    double orient2d(double*, double*, double*);
    double orient2dfast(double*, double*, double*);

#ifdef __cplusplus
}
#endif
