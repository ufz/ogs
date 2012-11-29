/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshQualityChecker.h
 *
 *  Created on 2010-12-08 by Thomas Fischer
 */

#ifndef MESHQUALITYCHECKER_H_
#define MESHQUALITYCHECKER_H_

#include <vector>

// BaseLib
#include "Histogram.h"

// MSH
#include "Mesh.h"
#include "Elements/Element.h"

#include "logog/include/logog.hpp"

namespace MeshLib
{
class MeshQualityChecker
{
public:
	MeshQualityChecker(Mesh const* const mesh);

	virtual ~MeshQualityChecker () {}

	virtual void check () = 0;
	std::vector<double> const& getMeshQuality () const;
	double getMinValue() const;
	double getMaxValue() const;
	virtual BaseLib::Histogram<double> getHistogram (std::size_t nclasses = 0) const;

protected:
	void errorMsg (const Element* elem, std::size_t idx) const;

	double _min;
	double _max;
	Mesh const* const _mesh;
	std::vector<double> _mesh_quality_measure;
};
}

#endif /* MESHQUALITYCHECKER_H_ */
