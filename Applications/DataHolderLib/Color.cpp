/**
 * \file
 * \author Karsten Rink
 * \date   2010-06-16
 * \brief  Implementation of the Color class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Color.h"

#include "BaseLib/Logging.h"

namespace DataHolderLib
{
Color createColor(unsigned char r,
                  unsigned char g,
                  unsigned char b,
                  unsigned char a)
{
    return Color{{r, g, b, a}};
}

Color getRandomColor()
{
    return createColor(static_cast<unsigned char>((rand() % 5) * 50),
                       static_cast<unsigned char>((rand() % 5) * 50),
                       static_cast<unsigned char>((rand() % 5) * 50));
}

Color getColor(const std::string& id, std::map<std::string, Color>& colors)
{
    auto it = colors.find(id);

    if (it == end(colors))
    {
        WARN("Key '{:s}' not found in color lookup table.", id);
        it = colors.insert({id, getRandomColor()}).first;
    }

    return it->second;
}

}  // namespace DataHolderLib
