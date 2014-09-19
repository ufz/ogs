/**
 * @file
 * @date 2014-09-19
 * @brief Implementation of heuristic search length strategy.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "HeuristicSearchLength.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"

namespace MeshGeoToolsLib
{

HeuristicSearchLength::HeuristicSearchLength(MeshLib::Mesh const& mesh)
: SearchLength(mesh)
{
	double sum (0.0);
	double sum_of_sqr (0.0);
	std::size_t edge_cnt(0);
	std::vector<MeshLib::Element*> const& elements(_mesh.getElements());

	for (std::vector<MeshLib::Element*>::const_iterator it(elements.cbegin());
			it != elements.cend(); ++it) {
		std::size_t const n_edges((*it)->getNEdges());
		for (std::size_t k(0); k<n_edges; k++) {
			MeshLib::Line const* edge(static_cast<MeshLib::Line const*>((*it)->getEdge(k)));
			if (!edge) {
				delete edge;
				continue;
			}
			double const len(edge->getLength());
			sum += len;
			sum_of_sqr += len*len;
			delete edge;
		}
		edge_cnt += n_edges;
	}

	const double mu (sum/edge_cnt);
	const double s (sqrt(1.0/(edge_cnt-1) * (sum_of_sqr - (sum*sum)/edge_cnt) ));
	// heuristic to prevent negative search lengths
	// in the case of a big standard deviation s
	double c(2.0);
	while (mu < c * s) {
		c *= 0.9;
	}

	_search_length = (mu - c * s)/2;

	DBUG("[MeshNodeSearcher::MeshNodeSearcher] Calculated search length for mesh \"%s\" is %f.",
			_mesh.getName().c_str(), _search_length);
}

double HeuristicSearchLength::getSearchLength() const
{
	return _search_length;
}

} // end namespace MeshGeoToolsLib
