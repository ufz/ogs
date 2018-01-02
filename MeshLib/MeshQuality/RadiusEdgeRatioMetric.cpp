/**
 * \file   RadiusEdgeRatioMetric.cpp
 * \author Karsten Rink
 * \date   2014-09-02
 * \brief  Implementation of the RadiusEdgeRadioMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RadiusEdgeRatioMetric.h"

#include "GeoLib/MinimalBoundingSphere.h"
#include "MeshLib/Node.h"

namespace MeshLib
{

RadiusEdgeRatioMetric::RadiusEdgeRatioMetric(Mesh const& mesh)
: ElementQualityMetric(mesh)
{}

void RadiusEdgeRatioMetric::calculateQuality ()
{
    std::vector<MeshLib::Element*> const& elements(_mesh.getElements());
    std::size_t const nElements (_mesh.getNumberOfElements());
    for (std::size_t k(0); k < nElements; k++)
    {
        Element const& elem (*elements[k]);
        std::size_t const n_nodes (elem.getNumberOfBaseNodes());
        std::vector<MathLib::Point3d*> pnts(n_nodes);
        std::copy_n(elem.getNodes(), n_nodes, pnts.begin());
        GeoLib::MinimalBoundingSphere const s(pnts);
        double min, max;
        elem.computeSqrEdgeLengthRange(min, max);
        _element_quality_metric[k] = sqrt(min)/(2*s.getRadius());
    }
}

} // end namespace MeshLib
