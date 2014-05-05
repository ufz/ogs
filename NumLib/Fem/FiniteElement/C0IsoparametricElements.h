/**
 * \author Norihiro Watanabe
 * \date   2013-08-13
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef C0ISOPARAMETRICELEMENTS_H_
#define C0ISOPARAMETRICELEMENTS_H_

#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"
#include "NumLib/Fem/Integration/IntegrationGaussTri.h"

#include "TemplateIsoparametric.h"

namespace NumLib
{

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct FeLINE2
{
    typedef TemplateIsoparametric<MeshLib::Line, ShapeLine2, IntegrationGaussRegular<1>, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct FeTRI3
{
    typedef TemplateIsoparametric<MeshLib::Tri, ShapeTri3, IntegrationGaussTri, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct FeQUAD4
{
    typedef TemplateIsoparametric<MeshLib::Quad, ShapeQuad4, IntegrationGaussRegular<2>, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct FeHEX8
{
    typedef TemplateIsoparametric<MeshLib::Hex, ShapeHex8, IntegrationGaussRegular<3>, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

} // NumLib

#endif //C0ISOPARAMETRICELEMENTS_H_
