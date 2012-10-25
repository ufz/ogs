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
    int x = 0x00000001;
    if (*(char*)&x)
        return true;              //am little
    else
        return false;             //am big
}

}

#endif
