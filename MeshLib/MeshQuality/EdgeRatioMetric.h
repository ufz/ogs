/**
 * \file   EdgeRatioMetric.h
 * \author Thomas Fischer
 * \date   2011-03-03
 * \brief  Definition of the AreaMetric class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGERATIOMETRIC_H_
#define EDGERATIOMETRIC_H_

#include "ElementQualityMetric.h"
#include "GeoLib/Point.h"

namespace MeshLib
{

/** 
 * Calculates the quality of mesh elements based on the ratio between shortest and longest edge of an element
 */
class EdgeRatioMetric : public ElementQualityMetric
{
public:
	EdgeRatioMetric(Mesh const& mesh);
	virtual ~EdgeRatioMetric () {}

	virtual void calculateQuality ();

private:
	double checkTriangle (GeoLib::Point const& a,
	                      GeoLib::Point const& b,
	                      GeoLib::Point const& c) const;
	double checkQuad (GeoLib::Point const& a,
	                  GeoLib::Point const& b,
	                  GeoLib::Point const& c,
	                  GeoLib::Point const& d) const;
	double checkTetrahedron (GeoLib::Point const& a,
	                         GeoLib::Point const& b,
	                         GeoLib::Point const& c,
	                         GeoLib::Point const& d) const;
	double checkPrism (std::vector<const GeoLib::Point*> const& pnts) const;
	double checkPyramid (std::vector<const GeoLib::Point*> const& pnts) const;
	double checkHexahedron (std::vector<const GeoLib::Point*> const& pnts) const;
};
}

#endif /* EDGERATIOMETRIC_H_ */
