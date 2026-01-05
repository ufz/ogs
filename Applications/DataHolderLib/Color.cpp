// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
