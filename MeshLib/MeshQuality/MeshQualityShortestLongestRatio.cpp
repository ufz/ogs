/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshQualityShortestLongestRatio.cpp
 *
 *  Created on 2011-03-03 by Thomas Fischer
 */

#include "MeshQualityShortestLongestRatio.h"
#include "Node.h"
#include "MathTools.h"

namespace MeshLib
{
MeshQualityShortestLongestRatio::MeshQualityShortestLongestRatio(
        Mesh const* const mesh) :
	MeshQualityChecker(mesh)
{
}

void MeshQualityShortestLongestRatio::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::Element*>& elements(_mesh->getElements());
	const size_t nElements (_mesh->getNElements());
	for (size_t k(0); k < nElements; k++)
	{
		const Element* elem (elements[k]);
		switch (elem->getGeomType())
		{
		case MshElemType::EDGE:
			_mesh_quality_measure[k] = 1.0;
			break;
		case MshElemType::TRIANGLE: {
			_mesh_quality_measure[k] = checkTriangle(elem->getNode(0), elem->getNode(1), elem->getNode(2));
			break;
		}
		case MshElemType::QUAD: {
			_mesh_quality_measure[k] = checkQuad(elem->getNode(0), elem->getNode(1), elem->getNode(2), elem->getNode(3));
			break;
		}
		case MshElemType::TETRAHEDRON: {
			_mesh_quality_measure[k] = checkTetrahedron(elem->getNode(0), elem->getNode(1), elem->getNode(2), elem->getNode(3));
			break;
		}
		case MshElemType::PRISM: {
			std::vector<const GeoLib::Point*> pnts;
			for (size_t j(0); j < 6; j++)
				pnts.push_back(elem->getNode(j));
			_mesh_quality_measure[k] = checkPrism(pnts);
			break;
		}
		case MshElemType::HEXAHEDRON: {
			std::vector<const GeoLib::Point*> pnts;
			for (size_t j(0); j < 8; j++)
				pnts.push_back(elem->getNode(j));
			_mesh_quality_measure[k] = checkHexahedron(pnts);
			break;
		}
		default:
			ERR ("MeshQualityShortestLongestRatio::check () check for element type %s not implemented.",
			     MshElemType2String(elem->getGeomType()).c_str());
		}
	}
}

double MeshQualityShortestLongestRatio::checkTriangle (GeoLib::Point const* const a,
                                                       GeoLib::Point const* const b,
                                                       GeoLib::Point const* const c) const
{
	double len0 (sqrt(MathLib::sqrDist (b,a)));
	double len1 (sqrt(MathLib::sqrDist (b,c)));
	double len2 (sqrt(MathLib::sqrDist (a,c)));

	if (len0 < len1 && len0 < len2)
	{
		if (len1 < len2)
			return len0 / len2;
		else
			return len0 / len1;
	}
	else
	{
		if (len1 < len2)
		{
			if (len0 < len2)
				return len1 / len2;
			else
				return len1 / len0;
		}
		else
		{
			if (len0 < len1)
				return len2 / len1;
			else
				return len2 / len0;
		}
	}
}

double MeshQualityShortestLongestRatio::checkQuad (GeoLib::Point const* const a,
                                                   GeoLib::Point const* const b,
                                                   GeoLib::Point const* const c,
                                                   GeoLib::Point const* const d) const
{
	double sqr_lengths[4] = {MathLib::sqrDist (b,a),
		                 MathLib::sqrDist (c,b),
		                 MathLib::sqrDist (d,c),
		                 MathLib::sqrDist (a,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 4; i++)
		for (size_t j(i + 1); j < 4; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[3]);
}

double MeshQualityShortestLongestRatio::checkTetrahedron (GeoLib::Point const* const a,
                                                          GeoLib::Point const* const b,
                                                          GeoLib::Point const* const c,
                                                          GeoLib::Point const* const d) const
{
	double sqr_lengths[6] = {MathLib::sqrDist (b,a), MathLib::sqrDist (c,b),
		                 MathLib::sqrDist (c,a), MathLib::sqrDist (a,d),
		                 MathLib::sqrDist (b,d), MathLib::sqrDist (c,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 6; i++)
		for (size_t j(i + 1); j < 6; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[5]);
}

double MeshQualityShortestLongestRatio::checkPrism (std::vector<const GeoLib::Point*> const & pnts) const
{
	double sqr_lengths[9] = {MathLib::sqrDist (pnts[0],pnts[1]),
		                 MathLib::sqrDist (pnts[1],pnts[2]),
		                 MathLib::sqrDist (pnts[2],pnts[0]),
		                 MathLib::sqrDist (pnts[3],pnts[4]),
		                 MathLib::sqrDist (pnts[4],pnts[5]),
		                 MathLib::sqrDist (pnts[5],pnts[3]),
		                 MathLib::sqrDist (pnts[0],pnts[3]),
		                 MathLib::sqrDist (pnts[1],pnts[4]),
		                 MathLib::sqrDist (pnts[2],pnts[5])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 9; i++)
		for (size_t j(i + 1); j < 9; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[8]);
}

double MeshQualityShortestLongestRatio::checkHexahedron (std::vector<const GeoLib::Point*> const & pnts)
const
{
	double sqr_lengths[12] = {MathLib::sqrDist (pnts[0],pnts[1]),
		                  MathLib::sqrDist (pnts[1],pnts[2]),
		                  MathLib::sqrDist (pnts[2],pnts[3]),
		                  MathLib::sqrDist (pnts[3],pnts[0]),
		                  MathLib::sqrDist (pnts[4],pnts[5]),
		                  MathLib::sqrDist (pnts[5],pnts[6]),
		                  MathLib::sqrDist (pnts[6],pnts[7]),
		                  MathLib::sqrDist (pnts[7],pnts[4]),
		                  MathLib::sqrDist (pnts[0],pnts[4]),
		                  MathLib::sqrDist (pnts[1],pnts[5]),
		                  MathLib::sqrDist (pnts[2],pnts[6]),
		                  MathLib::sqrDist (pnts[3],pnts[7])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i < 12; i++)
		for (size_t j(i + 1); j < 12; j++)
			if (sqr_lengths[i] >= sqr_lengths[j])
				std::swap (sqr_lengths[i], sqr_lengths[j]);

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[11]);
}
} // end namespace MeshLib
