/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NaturalNodeCoordinates.h"

namespace NumLib
{
constexpr std::array<std::array<double, 3>, 2>
    NaturalCoordinates<MeshLib::Line>::coordinates;

constexpr std::array<std::array<double, 3>, 3>
    NaturalCoordinates<MeshLib::Line3>::coordinates;

constexpr std::array<std::array<double, 3>, 3>
    NaturalCoordinates<MeshLib::Tri>::coordinates;

constexpr std::array<std::array<double, 3>, 6>
    NaturalCoordinates<MeshLib::Tri6>::coordinates;

constexpr std::array<std::array<double, 3>, 4>
    NaturalCoordinates<MeshLib::Quad>::coordinates;

constexpr std::array<std::array<double, 3>, 8>
    NaturalCoordinates<MeshLib::Quad8>::coordinates;

constexpr std::array<std::array<double, 3>, 9>
    NaturalCoordinates<MeshLib::Quad9>::coordinates;

constexpr std::array<std::array<double, 3>, 4>
    NaturalCoordinates<MeshLib::Tet>::coordinates;

constexpr std::array<std::array<double, 3>, 10>
    NaturalCoordinates<MeshLib::Tet10>::coordinates;

constexpr std::array<std::array<double, 3>, 6>
    NaturalCoordinates<MeshLib::Prism>::coordinates;

constexpr std::array<std::array<double, 3>, 15>
    NaturalCoordinates<MeshLib::Prism15>::coordinates;

constexpr std::array<std::array<double, 3>, 5>
    NaturalCoordinates<MeshLib::Pyramid>::coordinates;

constexpr std::array<std::array<double, 3>, 13>
    NaturalCoordinates<MeshLib::Pyramid13>::coordinates;

constexpr std::array<std::array<double, 3>, 8>
    NaturalCoordinates<MeshLib::Hex>::coordinates;

constexpr std::array<std::array<double, 3>, 20>
    NaturalCoordinates<MeshLib::Hex20>::coordinates;
}  // namespace NumLib
