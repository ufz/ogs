/*
 * MeshQualityArea.cpp
 *
 * 2011/03/17 KR Initial Implementation
 */

#include "MeshQualityArea.h"
#include "MathTools.h"

namespace MeshLib
{
MeshQualityArea::MeshQualityArea(Mesh const* const mesh)
	: MeshQualityChecker(mesh)
{}

void MeshQualityArea::check()
{
	// get all elements of mesh
	const std::vector<MeshLib::Element*> &elements(_mesh->getElements());

	const size_t nElems(elements.size());
	for (size_t k(0); k < nElems; k++) 
	{
		double area(std::numeric_limits<double>::max());
		const Element* elem (elements[k]);

		if (elem->getDimension() == 1)
		{
			_mesh_quality_measure[k] = -1.0;
			continue;
		}
		else if (elem->getDimension() == 2)
		{		
			area = elem->getContent();
			if (area < sqrt(fabs(std::numeric_limits<double>::min()))) errorMsg(elem, k);
		} 
		else {
			size_t nFaces(elem->getNFaces());

			for (size_t i = 0; i < nFaces; i++) 
			{
				const double sub_area (elem->getFace(i)->getContent());

				if (sub_area < sqrt(fabs(std::numeric_limits<double>::min())))
					errorMsg(elem, k);
				if (sub_area < area) area = sub_area;
			}
		}
		// update _min and _max values
		if (_min > area) _min = area;
		if (_max < area) _max = area;
		_mesh_quality_measure[k] = area;
	}
}

} // end namespace MeshLib
