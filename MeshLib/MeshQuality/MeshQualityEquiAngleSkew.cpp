/**
 * \file
 * \author Thomas Fischer
 * \date   2011-03-17
 * \brief  Implementation of the MeshQualityEquiAngleSkew class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshQualityEquiAngleSkew.h"
#include "Node.h"

#include "MathTools.h"

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

namespace MeshLib
{
MeshQualityEquiAngleSkew::MeshQualityEquiAngleSkew(Mesh const* const mesh) :
	MeshQualityChecker(mesh), M_PI_THIRD (M_PI / 3.0), TWICE_M_PI (2 * M_PI)
{}

MeshQualityEquiAngleSkew::~MeshQualityEquiAngleSkew()
{}

void MeshQualityEquiAngleSkew::check ()
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
			_mesh_quality_measure[k] = -1.0;
			break;
		case MshElemType::TRIANGLE:
			_mesh_quality_measure[k] = checkTriangle (elem);
			break;
		case MshElemType::QUAD:
			_mesh_quality_measure[k] = checkQuad (elem);
			break;
		case MshElemType::TETRAHEDRON:
			_mesh_quality_measure[k] = checkTetrahedron (elem);
			break;
		case MshElemType::HEXAHEDRON:
			_mesh_quality_measure[k] = checkHexahedron (elem);
			break;
		case MshElemType::PRISM:
			_mesh_quality_measure[k] = checkPrism (elem);
			break;
		default:
			break;
		}
	}
}

double MeshQualityEquiAngleSkew::checkTriangle (Element const* const elem) const
{
	double const* const node0 (elem->getNode(0)->getCoords());
	double const* const node1 (elem->getNode(1)->getCoords());
	double const* const node2 (elem->getNode(2)->getCoords());

	double min_angle (M_PI_2), max_angle (0.0);
	getMinMaxAngleFromTriangle (node0, node1, node2, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_THIRD) / (M_PI - M_PI_THIRD),
	                (M_PI_THIRD - min_angle) / (M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkQuad (Element const* const elem) const
{
	double const* const node0 (elem->getNode(0)->getCoords());
	double const* const node1 (elem->getNode(1)->getCoords());
	double const* const node2 (elem->getNode(2)->getCoords());
	double const* const node3 (elem->getNode(3)->getCoords());

	double min_angle (TWICE_M_PI);
	double max_angle (0.0);

	getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_2) / (M_PI - M_PI_2), (M_PI_2 - min_angle) / (M_PI_2));
}

double MeshQualityEquiAngleSkew::checkTetrahedron (Element const* const elem) const
{
	double const* const node0 (elem->getNode(0)->getCoords());
	double const* const node1 (elem->getNode(1)->getCoords());
	double const* const node2 (elem->getNode(2)->getCoords());
	double const* const node3 (elem->getNode(3)->getCoords());

	double min_angle (M_PI_2);
	double max_angle (0.0);

	// first triangle (0,1,2)
	getMinMaxAngleFromTriangle(node0, node1, node2, min_angle, max_angle);
	// second triangle (0,1,3)
	getMinMaxAngleFromTriangle(node0, node1, node3, min_angle, max_angle);
	// third triangle (0,2,3)
	getMinMaxAngleFromTriangle(node0, node2, node3, min_angle, max_angle);
	// fourth triangle (1,2,3)
	getMinMaxAngleFromTriangle(node1, node2, node3, min_angle, max_angle);

	return 1.0 - std::max((max_angle - M_PI_2) / (M_PI - M_PI_THIRD),
	                      (M_PI_THIRD - min_angle) / (M_PI_THIRD));
}

double MeshQualityEquiAngleSkew::checkHexahedron (Element const* const elem) const
{
	double const* const node0 (elem->getNode(0)->getCoords());
	double const* const node1 (elem->getNode(1)->getCoords());
	double const* const node2 (elem->getNode(2)->getCoords());
	double const* const node3 (elem->getNode(3)->getCoords());
	double const* const node4 (elem->getNode(4)->getCoords());
	double const* const node5 (elem->getNode(5)->getCoords());
	double const* const node6 (elem->getNode(6)->getCoords());
	double const* const node7 (elem->getNode(7)->getCoords());

	double min_angle (2 * M_PI);
	double max_angle (0.0);

	// first surface (0,1,2,3)
	getMinMaxAngleFromQuad (node0, node1, node2, node3, min_angle, max_angle);
	// second surface (0,3,7,4)
	getMinMaxAngleFromQuad (node0, node3, node7, node4, min_angle, max_angle);
	// third surface (4,5,6,7)
	getMinMaxAngleFromQuad (node4, node5, node6, node7, min_angle, max_angle);
	// fourth surface (5,1,2,6)
	getMinMaxAngleFromQuad (node5, node1, node2, node6, min_angle, max_angle);
	// fifth surface (5,1,0,4)
	getMinMaxAngleFromQuad (node5, node1, node0, node4, min_angle, max_angle);
	// sixth surface (6,2,3,7)
	getMinMaxAngleFromQuad (node6, node2, node3, node7, min_angle, max_angle);

	return 1.0 -
	       std::max((max_angle - M_PI_2) / (M_PI - M_PI_2), (M_PI_2 - min_angle) / (M_PI_2));
}

double MeshQualityEquiAngleSkew::checkPrism (Element const* const elem) const
{
	double const* const node0 (elem->getNode(0)->getCoords());
	double const* const node1 (elem->getNode(1)->getCoords());
	double const* const node2 (elem->getNode(2)->getCoords());
	double const* const node3 (elem->getNode(3)->getCoords());
	double const* const node4 (elem->getNode(4)->getCoords());
	double const* const node5 (elem->getNode(5)->getCoords());

	double min_angle_tri (2 * M_PI);
	double max_angle_tri (0.0);

	// first triangle (0,1,2)
	getMinMaxAngleFromTriangle (node0, node1, node2, min_angle_tri, max_angle_tri);
	// second surface (3,4,5)
	getMinMaxAngleFromTriangle (node3, node4, node5, min_angle_tri, max_angle_tri);

	double tri_criterion (1.0 - std::max((max_angle_tri - M_PI_2) / (M_PI - M_PI_THIRD),
	                                     (M_PI_THIRD - min_angle_tri) / (M_PI_THIRD)));

	double min_angle_quad (2 * M_PI);
	double max_angle_quad (0.0);
	// surface (0,3,4,1)
	getMinMaxAngleFromQuad (node0, node3, node4, node1, min_angle_quad, max_angle_quad);
	// surface (2,5,3,0)
	getMinMaxAngleFromQuad (node2, node5, node3, node0, min_angle_quad, max_angle_quad);
	// surface (1,2,5,4)
	getMinMaxAngleFromQuad (node1, node2, node5, node4, min_angle_quad, max_angle_quad);

	double quad_criterion (1.0 - std::max((max_angle_quad - M_PI_2) / (M_PI - M_PI_2),
	                                      (M_PI_2 - min_angle_quad) / (M_PI_2)));

	return std::min (tri_criterion, quad_criterion);
}

void MeshQualityEquiAngleSkew::getMinMaxAngleFromQuad (
        double const* const n0, double const* const n1,
        double const* const n2, double const* const n3,
        double &min_angle, double &max_angle) const
{
	double angle (MathLib::getAngle (n3, n0, n1));
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n0, n1, n2);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n1, n2, n3);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n2, n3, n0);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;
}

void MeshQualityEquiAngleSkew::getMinMaxAngleFromTriangle(double const* const n0,
                                                          double const* const n1,
                                                          double const* const n2,
                                                          double &min_angle,
                                                          double &max_angle) const
{
	double angle (MathLib::getAngle (n2, n0, n1));
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n0, n1, n2);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;

	angle = MathLib::getAngle (n1, n2, n0);
	if (angle < min_angle)
		min_angle = angle;
	if (angle > max_angle)
		max_angle = angle;
}
} // end namespace MeshLib
