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
	switch (type)
	{
		case VTK_LINE:
			return CellType::LINE2;
		case VTK_QUADRATIC_EDGE:
			return CellType::LINE3;
		case VTK_TRIANGLE:
			return CellType::TRI3;
		case VTK_QUADRATIC_TRIANGLE:
			return CellType::TRI6;
		case VTK_QUAD:
			return CellType::QUAD4;
		case VTK_QUADRATIC_QUAD:
			return CellType::QUAD8;
		case VTK_BIQUADRATIC_QUAD:
			return CellType::QUAD9;
		case VTK_HEXAHEDRON:
			return CellType::HEX8;
		case VTK_QUADRATIC_HEXAHEDRON:
			return CellType::HEX20;
		case VTK_TRIQUADRATIC_HEXAHEDRON:
			return CellType::HEX27;
		case VTK_TETRA:
			return CellType::TET4;
		case VTK_QUADRATIC_TETRA:
			return CellType::TET10;
		case VTK_WEDGE:
			return CellType::PRISM6;
		case VTK_QUADRATIC_WEDGE:
			return CellType::PRISM15;
		case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
			return CellType::PRISM18;
		case VTK_PYRAMID:
			return CellType::PYRAMID5;
		case VTK_QUADRATIC_PYRAMID:
			return CellType::PYRAMID13;
		default:
			return CellType::INVALID;
	}
}

int OGSToVtkCellType(CellType ogs)
{
	switch (ogs)
	{
		case CellType::LINE2:
			return VTK_LINE;
		case CellType::LINE3:
			return VTK_QUADRATIC_EDGE;
		case CellType::TRI3:
			return VTK_TRIANGLE;
		case CellType::TRI6:
			return VTK_QUADRATIC_TRIANGLE;
		case CellType::QUAD4:
			return VTK_QUAD;
		case CellType::QUAD8:
			return VTK_QUADRATIC_QUAD;
		case CellType::QUAD9:
			return VTK_BIQUADRATIC_QUAD;
		case CellType::HEX8:
			return VTK_HEXAHEDRON;
		case CellType::HEX20:
			return VTK_QUADRATIC_HEXAHEDRON;
		case CellType::HEX27:
			return VTK_TRIQUADRATIC_HEXAHEDRON;
		case CellType::TET4:
			return VTK_TETRA;
		case CellType::TET10:
			return VTK_QUADRATIC_TETRA;
		case CellType::PRISM6:
			return VTK_WEDGE;
		case CellType::PRISM15:
			return VTK_QUADRATIC_WEDGE;
		case CellType::PRISM18:
			return VTK_BIQUADRATIC_QUADRATIC_WEDGE;
		case CellType::PYRAMID5:
			return VTK_PYRAMID;
		case CellType::PYRAMID13:
			return VTK_QUADRATIC_PYRAMID;
		default:
			return -1;
	}
}


