/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "VtkOGSEnum.h"

#include <vtkCellType.h>

#include "BaseLib/Error.h"
#include "BaseLib/cpp23.h"

int OGSToVtkCellType(MeshLib::CellType ogs)
{
    switch (ogs)
    {
        case MeshLib::CellType::POINT1:
            return VTK_VERTEX;
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
            OGS_FATAL(
                "Unknown cell type in conversion from OGS to VTK. Given cell "
                "type value is {}.",
                BaseLib::to_underlying(ogs));
    }
}
