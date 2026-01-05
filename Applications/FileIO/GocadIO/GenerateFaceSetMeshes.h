// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>

#include "MeshLib/Mesh.h"

#include "GocadSGridReader.h"

namespace FileIO
{
namespace Gocad
{

void generateFaceSets(GocadSGridReader const& reader, std::string const& path);

}  //  namespace Gocad
}  //  namespace FileIO
