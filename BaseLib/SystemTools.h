/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file SystemTools.h
 *
 * Created on 2012-07-16 by Norihiro Watanabe
 */

#pragma once


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
