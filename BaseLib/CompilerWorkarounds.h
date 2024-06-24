/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#if defined(__GNUC__) && __GNUC__ == 14
#define OGS_NO_DANGLING [[gnu::no_dangling]]
#else
#define OGS_NO_DANGLING
#endif
