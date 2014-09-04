/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsSearcher.h"

#include "GeoLib/GeoObject.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSearcher.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsAlongPolyline.h"
#include "MeshGeoToolsLib/BoundaryElementsAlongSurface.h"


namespace MeshGeoToolsLib
{

BoundaryElementsSearcher::BoundaryElementsSearcher(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher) : _mesh(mesh), _mshNodeSearcher(mshNodeSearcher)
{}

BoundaryElementsSearcher::~BoundaryElementsSearcher()
{
	for (auto p : _boundary_elements_along_polylines)
		delete p;
	for (auto p : _boundary_elements_along_surfaces)
		delete p;
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElements(GeoLib::GeoObject const& geoObj)
{
	switch (geoObj.getGeoType()) {
	case GeoLib::GEOTYPE::POLYLINE:
		return this->getBoundaryElementsAlongPolyline(*dynamic_cast<const GeoLib::Polyline*>(&geoObj));
		break;
	case GeoLib::GEOTYPE::SURFACE:
		return this->getBoundaryElementsAlongSurface(*dynamic_cast<const GeoLib::Surface*>(&geoObj));
		break;
	default:
		const static std::vector<MeshLib::Element*> dummy;
		return dummy;
	}
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElementsAlongPolyline(GeoLib::Polyline const& ply)
{
	std::vector<BoundaryElementsAlongPolyline*>::const_iterator it(_boundary_elements_along_polylines.begin());
	for (; it != _boundary_elements_along_polylines.end(); ++it) {
		if (&(*it)->getPolyline() == &ply) {
			// we calculated mesh nodes for this polyline already
			return (*it)->getBoundaryElements();
		}
	}

	_boundary_elements_along_polylines.push_back(
			new BoundaryElementsAlongPolyline(_mesh, _mshNodeSearcher, ply));
	return _boundary_elements_along_polylines.back()->getBoundaryElements();
}

std::vector<MeshLib::Element*> const& BoundaryElementsSearcher::getBoundaryElementsAlongSurface(GeoLib::Surface const& sfc)
{
	std::vector<BoundaryElementsAlongSurface*>::const_iterator it(_boundary_elements_along_surfaces.begin());
	for (; it != _boundary_elements_along_surfaces.end(); ++it) {
		if (&(*it)->getSurface() == &sfc) {
			// we calculated mesh nodes for this surface already
			return (*it)->getBoundaryElements();
		}
	}

	_boundary_elements_along_surfaces.push_back(
			new BoundaryElementsAlongSurface(_mesh, _mshNodeSearcher, sfc));
	return _boundary_elements_along_surfaces.back()->getBoundaryElements();
}

} // end namespace MeshGeoTools

