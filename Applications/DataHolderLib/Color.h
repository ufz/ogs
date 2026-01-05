// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <array>
#include <map>
#include <string>

namespace DataHolderLib
{
using Color = std::array<unsigned char, 4>;

Color createColor(unsigned char r,
                  unsigned char g,
                  unsigned char b,
                  unsigned char a = 255);

/// Returns a random RGB colour.
Color getRandomColor();

/// Uses a color-lookup-table (in form of a map) to return a colour for a specified name. If the name is not
/// in the colortable a new entry is created with the new name and a random colour.
Color getColor(const std::string& id,
               std::map<std::string, DataHolderLib::Color>& colors);

}  // namespace DataHolderLib
