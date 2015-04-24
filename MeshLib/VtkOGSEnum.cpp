/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkOGSEnum.h"

#include <vtkCellType.h>

CellType VtkCellTypeToOGS(int type)
{
	CellType ogs;
	switch (type)
	{
		case VTK_LINE:
			ogs = CellType::LINE2;
			break;
		case VTK_QUADRATIC_EDGE:
			ogs = CellType::LINE3;
			break;
		case VTK_TRIANGLE:
			ogs = CellType::TRI3;
			break;
		case VTK_QUADRATIC_TRIANGLE:
			ogs = CellType::TRI6;
			break;
		case VTK_QUAD:
			ogs = CellType::QUAD4;
			break;
		case VTK_QUADRATIC_QUAD:
			ogs = CellType::QUAD8;
			break;
		case VTK_BIQUADRATIC_QUAD:
			ogs = CellType::QUAD9;
			break;
		case VTK_HEXAHEDRON:
			ogs = CellType::HEX8;
			break;
		case VTK_QUADRATIC_HEXAHEDRON:
			ogs = CellType::HEX20;
			break;
		case VTK_TRIQUADRATIC_HEXAHEDRON:
			ogs = CellType::HEX27;
			break;
		case VTK_TETRA:
			ogs = CellType::TET4;
			break;
		case VTK_QUADRATIC_TETRA:
			ogs = CellType::TET10;
			break;
		case VTK_WEDGE:
			ogs = CellType::PRISM6;
			break;
		case VTK_QUADRATIC_WEDGE:
			ogs = CellType::PRISM15;
			break;
		case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
			ogs = CellType::PRISM18;
			break;
		case VTK_PYRAMID:
			ogs = CellType::PYRAMID5;
			break;
		case VTK_QUADRATIC_PYRAMID:
			ogs = CellType::PYRAMID13;
			break;
		default:
			ogs = CellType::INVALID;
			break;
	}
	return ogs;
}

int OGSToVtkCellType(CellType ogs)
{
	int type = 0;
	switch (ogs)
	{
		case CellType::INVALID:
			break;
		case CellType::LINE2:
			type = VTK_LINE;
			break;
		case CellType::LINE3:
			type = VTK_QUADRATIC_EDGE;
			break;
		case CellType::TRI3:
			type = VTK_TRIANGLE;
			break;
		case CellType::TRI6:
			type = VTK_QUADRATIC_TRIANGLE;
			break;
		case CellType::QUAD4:
			type = VTK_QUAD;
			break;
		case CellType::QUAD8:
			type = VTK_QUADRATIC_QUAD;
			break;
		case CellType::QUAD9:
			type = VTK_BIQUADRATIC_QUAD;
			break;
		case CellType::HEX8:
			type = VTK_HEXAHEDRON;
			break;
		case CellType::HEX20:
			type = VTK_QUADRATIC_HEXAHEDRON;
			break;
		case CellType::HEX27:
			type = VTK_TRIQUADRATIC_HEXAHEDRON;
			break;
		case CellType::TET4:
			type = VTK_TETRA;
			break;
		case CellType::TET10:
			type = VTK_QUADRATIC_TETRA;
			break;
		case CellType::PRISM6:
			type = VTK_WEDGE;
			break;
		case CellType::PRISM15:
			type = VTK_QUADRATIC_WEDGE;
			break;
		case CellType::PRISM18:
			type = VTK_BIQUADRATIC_QUADRATIC_WEDGE;
			break;
		case CellType::PYRAMID5:
			type = VTK_PYRAMID;
			break;
		case CellType::PYRAMID13:
			type = VTK_QUADRATIC_PYRAMID;
			break;
	}
	return type;
}


