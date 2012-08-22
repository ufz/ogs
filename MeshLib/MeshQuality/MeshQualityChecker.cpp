/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "Node.h"
#include "Point.h"
#include <cmath>
#include <iostream>

namespace MeshLib
{
MeshQualityChecker::MeshQualityChecker(Mesh const* const mesh) :
	_min (-1.0), _max (-1.0), _mesh (mesh)
{
	if (_mesh)
		_mesh_quality_measure.resize (_mesh->getNElements(), -1.0);
}

BASELIB::Histogram<double> MeshQualityChecker::getHistogram (size_t nclasses) const
{
	if (nclasses == 0) {
		// simple suggestion: number of classes with Sturges criterion
		nclasses = static_cast<size_t>(1 + 3.3 * log (static_cast<float>((_mesh->getNElements()))));
	}

	return BASELIB::Histogram<double>(getMeshQuality(), nclasses, true);
}

void MeshQualityChecker::errorMsg (const Element* elem, size_t idx) const
{
	std::cout << "Error in MeshQualityChecker::check() - "
			  << "Calculated value of element is below double precision minimum." << std::endl;
	std::cout << "Points of " << MshElemType2String(elem->getType()) << "-Element " << idx << ": " << std::endl;
	for (size_t i(0); i < elem->getNNodes(); i++)
		std::cout << "\t Node " << i << " " << GeoLib::Point((elem->getNode(i))->getCoords()) << std::endl;
}

std::vector<double> const&
MeshQualityChecker::getMeshQuality () const
{
	return _mesh_quality_measure;
}

double MeshQualityChecker::getMinValue() const
{
	return _min;
}

double MeshQualityChecker::getMaxValue() const
{
	return _max;
}

} // end namespace MeshLib
