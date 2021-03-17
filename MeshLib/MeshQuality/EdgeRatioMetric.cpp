/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Implementation of the EdgeRatioMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EdgeRatioMetric.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
EdgeRatioMetric::EdgeRatioMetric(Mesh const& mesh) :
    ElementQualityMetric(mesh)
{
}

void EdgeRatioMetric::calculateQuality()
{
    auto const& elements(_mesh.getElements());
    auto const n_elements(_mesh.getNumberOfElements());
    for (std::size_t k(0); k < n_elements; k++)
    {
        auto const& elem(*elements[k]);
        std::unique_ptr<Element const> first_edge{elem.getEdge(0)};
        auto sqr_min_edge_length =
            MathLib::sqrDist(*first_edge->getNode(1), *first_edge->getNode(0));
        auto sqr_max_edge_length = sqr_min_edge_length;
        auto const n_edges(elem.getNumberOfEdges());
        for (std::size_t i = 1; i < n_edges; i++)
        {
            std::unique_ptr<Element const> edge{elem.getEdge(i)};
            auto const sqr_edge_length =
                MathLib::sqrDist(*edge->getNode(1), *edge->getNode(0));
            if (sqr_edge_length < sqr_min_edge_length)
            {
                sqr_min_edge_length = sqr_edge_length;
            }
            if (sqr_edge_length > sqr_max_edge_length)
            {
                sqr_max_edge_length = sqr_edge_length;
            }
        }
        _element_quality_metric[k] =
            std::sqrt(sqr_min_edge_length / sqr_max_edge_length);
    }
}

} // end namespace MeshLib
