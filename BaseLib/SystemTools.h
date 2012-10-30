/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file SystemTools.h
 *
 * Created on 2012-07-16 by Norihiro Watanabe
 */

#ifndef SYSTEMTOOLS_H
#define SYSTEMTOOLS_H

namespace BaseLib
{

/// return if this system supports little endian or not
inline bool IsLittleEndian()
{
#ifdef ENDIAN_IS_BIG
    return false;
#elif ENDIAN_IS_LITTLE
    return true;
#endif
}

}

#endif
