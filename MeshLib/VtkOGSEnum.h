/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKOGSENUM_H_
#define VTKOGSENUM_H_

#include "MeshEnums.h"

MeshLib::CellType VtkCellTypeToOGS(int type);

int OGSToVtkCellType(MeshLib::CellType ogs);

#endif // VTKOGSENUM_H_
