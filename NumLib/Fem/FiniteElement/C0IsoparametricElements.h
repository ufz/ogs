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

#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/Integration/IntegrationGaussRegular.h"

#include "TemplateIsoparametric.h"

namespace NumLib
{

template <class T_SHAPE_VECTOR, class T_DSHAPE_MATRIX, class T_JACOBIAN_MATRIX>
struct FeQUAD4
{
    typedef TemplateIsoparametric<MeshLib::Quad, ShapeQuad4, IntegrationGaussRegular<2>, T_SHAPE_VECTOR, T_DSHAPE_MATRIX, T_JACOBIAN_MATRIX> type;
};

} // NumLib

#endif //C0ISOPARAMETRICELEMENTS_H_
