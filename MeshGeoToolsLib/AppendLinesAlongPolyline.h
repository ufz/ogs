// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace GeoLib
{
class Polyline;
template <typename T> class TemplateVec;
using PolylineVec = TemplateVec<GeoLib::Polyline>;
}

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{

/**
 * Add line elements to a copy of a given mesh
 *
 * The function creates line elements from nodes located along user-provided
 * polylines. New elements will have a distinct material ID for each polyline.
 *
 * \remark The function allows creation of duplicated line elements.
 * \remark Line elements may not be placed along edges of existing elements.
 *
 * \param mesh      original mesh
 * \param ply_vec   polyline vector whose nodes are used to create line elements
 * \return a new mesh which is copied from a given mesh and additionally
 * includes line elements
 */
std::unique_ptr<MeshLib::Mesh> appendLinesAlongPolylines(
    const MeshLib::Mesh& mesh, const GeoLib::PolylineVec& ply_vec);
}  // namespace MeshGeoToolsLib
