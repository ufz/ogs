/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-08-13
 * \brief
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

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"
#include "NumLib/Fem/Extrapolation/ExtrapolationGaussQuad.h"

#include "TemplateIsoparametric.h"

namespace NumLib
{

#ifdef CXX_TEMPLATE_ALIASES_SUPPORTED

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
using QUAD4 = TemplateIsoparametric<MeshLib::Quad, ShapeQuad4, IntegrationGaussRegular<2>, ExtrapolationGaussQuad, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX>;

#else

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct QUAD4
{
    typedef TemplateIsoparametric<MeshLib::Quad, ShapeQuad4, IntegrationGaussRegular<2>, ExtrapolationGaussQuad, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

#endif

}

#endif //C0ISOPARAMETRICELEMENTS_H_
