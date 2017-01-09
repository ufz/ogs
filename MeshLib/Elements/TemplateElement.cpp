/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TemplateElement.h"

#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Pyramid.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Tet.h"
#include "MeshLib/Elements/Tri.h"

#ifndef WIN32
/// \todo Windows compiler does not accept this definition and issues a linking error.
template <class ELEMENT_RULE>
const unsigned MeshLib::TemplateElement<ELEMENT_RULE>::n_all_nodes;

template <class ELEMENT_RULE>
const unsigned MeshLib::TemplateElement<ELEMENT_RULE>::n_base_nodes;

template <class ELEMENT_RULE>
const unsigned MeshLib::TemplateElement<ELEMENT_RULE>::dimension;
#endif // WIN32

template class MeshLib::TemplateElement<MeshLib::HexRule20>;
template class MeshLib::TemplateElement<MeshLib::HexRule8>;
template class MeshLib::TemplateElement<MeshLib::LineRule2>;
template class MeshLib::TemplateElement<MeshLib::LineRule3>;
template class MeshLib::TemplateElement<MeshLib::PointRule1>;
template class MeshLib::TemplateElement<MeshLib::PrismRule15>;
template class MeshLib::TemplateElement<MeshLib::PrismRule6>;
template class MeshLib::TemplateElement<MeshLib::PyramidRule13>;
template class MeshLib::TemplateElement<MeshLib::PyramidRule5>;
template class MeshLib::TemplateElement<MeshLib::QuadRule4>;
template class MeshLib::TemplateElement<MeshLib::QuadRule8>;
template class MeshLib::TemplateElement<MeshLib::QuadRule9>;
template class MeshLib::TemplateElement<MeshLib::TetRule10>;
template class MeshLib::TemplateElement<MeshLib::TetRule4>;
template class MeshLib::TemplateElement<MeshLib::TriRule3>;
template class MeshLib::TemplateElement<MeshLib::TriRule6>;
