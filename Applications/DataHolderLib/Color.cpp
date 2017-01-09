/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-16
 * \brief  Implementation of the Color class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Color.h"

#include <fstream>
#include <sstream>

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"

namespace DataHolderLib {

Color createColor(unsigned char r, unsigned char g, unsigned char b)
{
    return Color{{r,g,b,255}};
}

Color createColor(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
{
    return Color{{r,g,b,a}};
}

Color getRandomColor()
{
    Color col;
    col[0] = static_cast<unsigned char>((rand()%5)*50);
    col[1] = static_cast<unsigned char>((rand()%5)*50);
    col[2] = static_cast<unsigned char>((rand()%5)*50);
    return col;
}

Color const getColor(const std::string &id, std::map<std::string, Color> &colors)
{
    for (std::map<std::string, Color>::const_iterator it=colors.begin(); it !=colors.end(); ++it)
    {
        if (id.compare(it->first) == 0)
            return it->second;
    }
    WARN("Key \"%s\" not found in color lookup table.", id.c_str());
    Color c = getRandomColor();
    colors.insert(std::pair<std::string, Color>(id, c));
    return c;
}

}
