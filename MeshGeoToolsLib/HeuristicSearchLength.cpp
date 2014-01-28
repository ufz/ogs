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
	double sum_of_sqr (0.0); // total length of edges
	std::size_t edge_cnt(0); // total length of edges squared
	std::vector<MeshLib::Element*> const& elements(_mesh.getElements());

	for (std::vector<MeshLib::Element*>::const_iterator it(elements.cbegin());
			it != elements.cend(); ++it) {
		std::size_t const n_edges((*it)->getNEdges());
		for (std::size_t k(0); k<n_edges; k++) {
			double const len =
				static_cast<MeshLib::Line const*>((*it)->getEdge(k))->getLength();
			sum += len;
			sum_of_sqr += len*len;
		}
		edge_cnt += n_edges;
	}

	const double mean (sum/edge_cnt);
	const double variance ((sum_of_sqr - (sum*sum)/edge_cnt)/(edge_cnt-1));

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

double HeuristicSearchLength::getSearchLength() const
{
	return _search_length;
}

} // end namespace MeshGeoToolsLib
