// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MeshLib
{
class Element;
}

namespace MeshToolsLib
{

double computeElementVolumeNumerically(MeshLib::Element const& e);
}  // namespace MeshToolsLib
