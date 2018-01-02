/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/ShapeFunction/ShapePoint1.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"

#include "TemplateIsoparametric.h"

namespace NumLib
{
template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FePOINT1
{
    using type =
        TemplateIsoparametric<ShapePoint1, T_SHAPE_MATRIX_POLICY<ShapePoint1>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeLINE2
{
    using type =
        TemplateIsoparametric<ShapeLine2, T_SHAPE_MATRIX_POLICY<ShapeLine2>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeLINE3
{
    using type =
        TemplateIsoparametric<ShapeLine3, T_SHAPE_MATRIX_POLICY<ShapeLine3>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeTRI3
{
    using type =
        TemplateIsoparametric<ShapeTri3, T_SHAPE_MATRIX_POLICY<ShapeTri3>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeTRI6
{
    using type =
        TemplateIsoparametric<ShapeTri6, T_SHAPE_MATRIX_POLICY<ShapeTri6>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeQUAD4
{
    using type =
        TemplateIsoparametric<ShapeQuad4, T_SHAPE_MATRIX_POLICY<ShapeQuad4>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeQUAD8
{
    using type =
        TemplateIsoparametric<ShapeQuad8, T_SHAPE_MATRIX_POLICY<ShapeQuad8>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeQUAD9
{
    using type =
        TemplateIsoparametric<ShapeQuad9, T_SHAPE_MATRIX_POLICY<ShapeQuad9>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeHEX8
{
    using type =
        TemplateIsoparametric<ShapeHex8, T_SHAPE_MATRIX_POLICY<ShapeHex8>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeHEX20
{
    using type =
        TemplateIsoparametric<ShapeHex20, T_SHAPE_MATRIX_POLICY<ShapeHex20>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeTET4
{
    using type =
        TemplateIsoparametric<ShapeTet4, T_SHAPE_MATRIX_POLICY<ShapeTet4>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FeTET10
{
    using type =
        TemplateIsoparametric<ShapeTet10, T_SHAPE_MATRIX_POLICY<ShapeTet10>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FePRISM6
{
    using type =
        TemplateIsoparametric<ShapePrism6, T_SHAPE_MATRIX_POLICY<ShapePrism6>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FePRISM15
{
    using type = TemplateIsoparametric<ShapePrism15,
                                       T_SHAPE_MATRIX_POLICY<ShapePrism15>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FePYRA5
{
    using type =
        TemplateIsoparametric<ShapePyra5, T_SHAPE_MATRIX_POLICY<ShapePyra5>>;
};

template <template <typename> class T_SHAPE_MATRIX_POLICY>
struct FePYRA13
{
    using type =
        TemplateIsoparametric<ShapePyra13, T_SHAPE_MATRIX_POLICY<ShapePyra13>>;
};

} // NumLib
