/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-09-22 13:28:09
 */

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif
#include <lis.h>
#ifdef conj
#undef conj
#endif
