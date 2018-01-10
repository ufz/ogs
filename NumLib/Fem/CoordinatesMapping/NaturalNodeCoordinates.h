/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cassert>

#include "BaseLib/Error.h"
#include "MathLib/Point3d.h"

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "MeshLib/Elements/HexRule20.h"
#include "MeshLib/Elements/HexRule8.h"
#include "MeshLib/Elements/LineRule2.h"
#include "MeshLib/Elements/LineRule3.h"
#include "MeshLib/Elements/PointRule1.h"
#include "MeshLib/Elements/PrismRule15.h"
#include "MeshLib/Elements/PrismRule6.h"
#include "MeshLib/Elements/PyramidRule13.h"
#include "MeshLib/Elements/PyramidRule5.h"
#include "MeshLib/Elements/QuadRule4.h"
#include "MeshLib/Elements/QuadRule8.h"
#include "MeshLib/Elements/QuadRule9.h"
#include "MeshLib/Elements/TemplateElement.h"
#include "MeshLib/Elements/TetRule10.h"
#include "MeshLib/Elements/TetRule4.h"
#include "MeshLib/Elements/TriRule3.h"
#include "MeshLib/Elements/TriRule6.h"

#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "ShapeMatrices.h"

namespace NumLib
{
/// A map 3D natural coordinates (-1 to 1, or 0 to 1) for each node of a
/// MeshLib::Element.
template <typename MeshLibElement>
struct NaturalCoordinates;

template <>
struct NaturalCoordinates<MeshLib::Line>
{
    static constexpr std::array<std::array<double, 3>, 2> coordinates = {
        {{{-1, 0, 0}}, {{1, 0, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Line3>
{
    static constexpr std::array<std::array<double, 3>, 3> coordinates = {
        {{{-1, 0, 0}}, {{1, 0, 0}}, {{0, 0, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Tri>
{
    static constexpr std::array<std::array<double, 3>, 3> coordinates = {
        {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Tri6>
{
    static constexpr std::array<std::array<double, 3>, 6> coordinates = {
        {{{0, 0, 0}},
         {{1, 0, 0}},
         {{0, 1, 0}},
         {{0.5, 0, 0}},
         {{0.5, 0.5, 0}},
         {{0, 0.5, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Quad>
{
    static constexpr std::array<std::array<double, 3>, 4> coordinates = {
        {{{1, 1, 0}}, {{-1, 1, 0}}, {{-1, -1, 0}}, {{1, -1, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Quad8>
{
    static constexpr std::array<std::array<double, 3>, 8> coordinates = {
        {{{1, 1, 0}},
         {{-1, 1, 0}},
         {{-1, -1, 0}},
         {{1, -1, 0}},
         {{0, 1, 0}},
         {{-1, 0, 0}},
         {{0, -1, 0}},
         {{1, 0, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Quad9>
{
    static constexpr std::array<std::array<double, 3>, 9> coordinates = {
        {{{1, 1, 0}},
         {{-1, 1, 0}},
         {{-1, -1, 0}},
         {{1, -1, 0}},
         {{0, 1, 0}},
         {{-1, 0, 0}},
         {{0, -1, 0}},
         {{1, 0, 0}},
         {{0, 0, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Tet>
{
    static constexpr std::array<std::array<double, 3>, 4> coordinates = {
        {{{0, 0, 0}}, {{1, 0, 0}}, {{0, 1, 0}}, {{0, 0, 1}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Tet10>
{
    static constexpr std::array<std::array<double, 3>, 10> coordinates = {
        {{{0, 0, 0}},
         {{1, 0, 0}},
         {{0, 1, 0}},
         {{0, 0, 1}},

         {{0.5, 0, 0}},
         {{0.5, 0.5, 0}},
         {{0, 0.5, 0}},

         {{0, 0, 0.5}},
         {{0.5, 0, 0.5}},
         {{0, 0.5, 0.5}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Prism>
{
    static constexpr std::array<std::array<double, 3>, 6> coordinates = {
        {{{0, 0, -1}},
         {{1, 0, -1}},
         {{0, 1, -1}},
         {{0, 0, 1}},
         {{1, 0, 1}},
         {{0, 1, 1}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Prism15>
{
    static constexpr std::array<std::array<double, 3>, 15> coordinates = {
        {{{0, 0, -1}},
         {{1, 0, -1}},
         {{0, 1, -1}},
         {{0, 0, 1}},
         {{1, 0, 1}},
         {{0, 1, 1}},

         {{0.5, 0, -1}},
         {{0.5, 0.5, -1}},
         {{0, 0.5, -1}},

         {{0.5, 0, 1}},
         {{0.5, 0.5, 1}},
         {{0, 0.5, 1}},

         {{0, 0, 0}},
         {{1, 0, 0}},
         {{0, 1, 0}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Pyramid>
{
    static constexpr std::array<std::array<double, 3>, 5> coordinates = {
        {{{-1, -1, -1}},
         {{1, -1, -1}},
         {{1, 1, -1}},
         {{-1, 1, -1}},
         {{0, 0, 1}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Pyramid13>
{
    static constexpr std::array<std::array<double, 3>, 13> coordinates = {{
        {{-1, -1, -1}},
        {{1, -1, -1}},
        {{1, 1, -1}},
        {{-1, 1, -1}},

        {{0, 0, 1}},

        {{0, -1, -1}},
        {{1, 0, -1}},
        {{0, 1, -1}},
        {{-1, 0, -1}},

        {{-1, -1, 0}},
        {{1, -1, 0}},
        {{1, 1, 0}},
        {{-1, 1, 0}},
    }};
};

template <>
struct NaturalCoordinates<MeshLib::Hex>
{
    static constexpr std::array<std::array<double, 3>,
                                MeshLib::Hex::n_all_nodes>
        coordinates = {{{{-1, -1, -1}},
                        {{1, -1, -1}},
                        {{1, 1, -1}},
                        {{-1, 1, -1}},
                        {{-1, -1, 1}},
                        {{1, -1, 1}},
                        {{1, 1, 1}},
                        {{-1, 1, 1}}}};
};

template <>
struct NaturalCoordinates<MeshLib::Hex20>
{
    static constexpr std::array<std::array<double, 3>,
                                MeshLib::Hex20::n_all_nodes>
        coordinates = {
            {{{-1, -1, -1}}, {{1, -1, -1}}, {{1, 1, -1}}, {{-1, 1, -1}},
             {{-1, -1, 1}},  {{1, -1, 1}},  {{1, 1, 1}},  {{-1, 1, 1}},

             {{0, -1, -1}},  {{1, 0, -1}},  {{0, 1, -1}}, {{-1, 0, -1}},
             {{0, -1, 1}},   {{1, 0, 1}},   {{0, 1, 1}},  {{-1, 0, 1}},

             {{-1, -1, 0}},  {{1, -1, 0}},  {{1, 1, 0}},  {{-1, 1, 0}}}};
};
}  // namespace NumLib
