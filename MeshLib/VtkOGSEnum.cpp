/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkOGSEnum.h"

#include <vtkCellType.h>

MeshLib::CellType VtkCellTypeToOGS(int type)
{
    switch (type)
    {
        case VTK_LINE:
            return MeshLib::CellType::LINE2;
        case VTK_QUADRATIC_EDGE:
            return MeshLib::CellType::LINE3;
        case VTK_TRIANGLE:
            return MeshLib::CellType::TRI3;
        case VTK_QUADRATIC_TRIANGLE:
            return MeshLib::CellType::TRI6;
        case VTK_QUAD:
            return MeshLib::CellType::QUAD4;
        case VTK_QUADRATIC_QUAD:
            return MeshLib::CellType::QUAD8;
        case VTK_BIQUADRATIC_QUAD:
            return MeshLib::CellType::QUAD9;
        case VTK_HEXAHEDRON:
            return MeshLib::CellType::HEX8;
        case VTK_QUADRATIC_HEXAHEDRON:
            return MeshLib::CellType::HEX20;
        case VTK_TRIQUADRATIC_HEXAHEDRON:
            return MeshLib::CellType::HEX27;
        case VTK_TETRA:
            return MeshLib::CellType::TET4;
        case VTK_QUADRATIC_TETRA:
            return MeshLib::CellType::TET10;
        case VTK_WEDGE:
            return MeshLib::CellType::PRISM6;
        case VTK_QUADRATIC_WEDGE:
            return MeshLib::CellType::PRISM15;
        case VTK_BIQUADRATIC_QUADRATIC_WEDGE:
            return MeshLib::CellType::PRISM18;
        case VTK_PYRAMID:
            return MeshLib::CellType::PYRAMID5;
        case VTK_QUADRATIC_PYRAMID:
            return MeshLib::CellType::PYRAMID13;
        default:
            return MeshLib::CellType::INVALID;
    }
}

int OGSToVtkCellType(MeshLib::CellType ogs)
{
    switch (ogs)
    {
        case MeshLib::CellType::LINE2:
            return VTK_LINE;
        case MeshLib::CellType::LINE3:
            return VTK_QUADRATIC_EDGE;
        case MeshLib::CellType::TRI3:
            return VTK_TRIANGLE;
        case MeshLib::CellType::TRI6:
            return VTK_QUADRATIC_TRIANGLE;
        case MeshLib::CellType::QUAD4:
            return VTK_QUAD;
        case MeshLib::CellType::QUAD8:
            return VTK_QUADRATIC_QUAD;
        case MeshLib::CellType::QUAD9:
            return VTK_BIQUADRATIC_QUAD;
        case MeshLib::CellType::HEX8:
            return VTK_HEXAHEDRON;
        case MeshLib::CellType::HEX20:
            return VTK_QUADRATIC_HEXAHEDRON;
        case MeshLib::CellType::HEX27:
            return VTK_TRIQUADRATIC_HEXAHEDRON;
        case MeshLib::CellType::TET4:
            return VTK_TETRA;
        case MeshLib::CellType::TET10:
            return VTK_QUADRATIC_TETRA;
        case MeshLib::CellType::PRISM6:
            return VTK_WEDGE;
        case MeshLib::CellType::PRISM15:
            return VTK_QUADRATIC_WEDGE;
        case MeshLib::CellType::PRISM18:
            return VTK_BIQUADRATIC_QUADRATIC_WEDGE;
        case MeshLib::CellType::PYRAMID5:
            return VTK_PYRAMID;
        case MeshLib::CellType::PYRAMID13:
            return VTK_QUADRATIC_PYRAMID;
        case MeshLib::CellType::INVALID:
            return -1;
        default:
            return -1;
    }
}


