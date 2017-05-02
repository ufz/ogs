/**
 * @date 2014-09-19
 * @brief Implementation of heuristic search length strategy.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "HeuristicSearchLength.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Elements/Element.h"

namespace MeshGeoToolsLib
{

HeuristicSearchLength::HeuristicSearchLength(MeshLib::Mesh const& mesh, LengthType length_type)
: _mesh(mesh)
{
    double sum (0.0);
    double sum_of_sqr (0.0); // total length of edges

    std::size_t n_sampling(0); // total length of edges squared
    std::vector<MeshLib::Element*> const& elements(_mesh.getElements());

    if (length_type==LengthType::Edge) {
        for (auto element : elements)
        {
            std::size_t const n_edges(element->getNumberOfEdges());
            for (std::size_t k(0); k<n_edges; k++) {
                auto edge =
                    element->getEdge(k);  // allocation inside getEdge().
                double const len = edge->getContent();
                delete edge;
                sum += len;
                sum_of_sqr += len*len;
            }
            n_sampling += n_edges;
        }
    } else {
        double min=0, max=0;
        for (const MeshLib::Element* e : elements) {
            e->computeSqrNodeDistanceRange(min, max, true);
            sum += std::sqrt(min);
            sum_of_sqr += min;
        }
        n_sampling = _mesh.getNumberOfElements();
    }

    const double mean (sum/n_sampling);
    const double variance ((sum_of_sqr - (sum*sum)/n_sampling)/(n_sampling-1));

    // Set the search length for the case of non-positive variance (which can
    // happen due to numerics).
    _search_length = mean/2;

    if (variance > 0) {
        if (variance < mean*mean/4)
            _search_length -= std::sqrt(variance);
        else
            _search_length = std::numeric_limits<double>::epsilon();
    }

    DBUG("[MeshNodeSearcher::MeshNodeSearcher] Calculated search length for mesh \"%s\" is %f.",
            _mesh.getName().c_str(), _search_length);
}

} // end namespace MeshGeoToolsLib
