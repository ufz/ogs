/*
 * MeshQualityChecker.h
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#ifndef MESHQUALITYCHECKER_H_
#define MESHQUALITYCHECKER_H_

#include <vector>

// BaseLib
#include "Histogram.h"

// MSH
#include "Mesh.h"
#include "Elements/Element.h"

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
	virtual BASELIB::Histogram<double> getHistogram (size_t nclasses = 0) const;

protected:
	void errorMsg (const Element* elem, size_t idx) const;

	double _min;
	double _max;
	Mesh const* const _mesh;
	std::vector<double> _mesh_quality_measure;
};
}

#endif /* MESHQUALITYCHECKER_H_ */
