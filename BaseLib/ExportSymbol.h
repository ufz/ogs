/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifndef OGS_EXPORT_SYMBOL
#if defined(WIN32) || defined(_WIN32)
#define OGS_EXPORT_SYMBOL __declspec(dllexport)
#else
#define OGS_EXPORT_SYMBOL __attribute__((visibility("default")))
#endif
#endif  // defined(OGS_EXPORT_SYMBOL)
