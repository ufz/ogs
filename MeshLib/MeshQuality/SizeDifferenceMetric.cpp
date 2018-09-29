/**
 * \file   SizeDifferenceMetric.cpp
 * \author Karsten Rink
 * \date   2015-03-24
 * \brief  Implementation of the SizeDifferenceMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SizeDifferenceMetric.h"

#include <limits>

namespace MeshLib
{

SizeDifferenceMetric::SizeDifferenceMetric(Mesh const& mesh) :
ElementQualityMetric(mesh)
{ }

void SizeDifferenceMetric::calculateQuality()
{
    std::vector<MeshLib::Element*> const& elements(_mesh.getElements());
    std::size_t const nElements (_mesh.getNumberOfElements());
    std::size_t const mesh_dim (_mesh.getDimension());

    for (std::size_t k=0; k < nElements; ++k)
    {
        Element const& elem (*elements[k]);
        if (elem.getDimension() < mesh_dim)
        {
            _element_quality_metric[k] = 0;
            continue;
        }

        std::size_t const n_neighbors (elem.getNumberOfNeighbors());
        double const vol_a (elem.getContent());

        double worst_ratio(1.0);
        for (std::size_t i=0; i < n_neighbors; ++i)
        {
            MeshLib::Element const*const neighbor (elem.getNeighbor(i));
            if (neighbor == nullptr)
                continue;
            double const vol_b (neighbor->getContent());
            double const ratio = (vol_a > vol_b) ? vol_b / vol_a : vol_a / vol_b;
            if (ratio < worst_ratio)
                worst_ratio = ratio;
        }
        _element_quality_metric[k] = worst_ratio;
    }
}

} // end namespace MeshLib
