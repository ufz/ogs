/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef BASELIB_ERROR_H
#define BASELIB_ERROR_H

#include <stdexcept>

#include "StringTools.h"
#include "FileTools.h"

#define OGS_STR(x) #x
#define OGS_STRINGIFY(x) OGS_STR(x)
#define OGS_LOCATION " at " + BaseLib::extractBaseName(__FILE__) + ", line " OGS_STRINGIFY(__LINE__)
#define OGS_FATAL(fmt, ...)\
    throw std::runtime_error(BaseLib::format(fmt, ##__VA_ARGS__) + OGS_LOCATION);

#endif //BASELIB_ERROR_H
