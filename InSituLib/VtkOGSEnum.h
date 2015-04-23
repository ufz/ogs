/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTKOGSENUM_H_
#define VTKOGSENUM_H_

#include "MeshEnums.h"

namespace InSituLib
{

CellType VtkCellTypeToOGS(int type);

int OGSToVtkCellType(CellType ogs);

} // end namespace

#endif // VTKOGSENUM_H_
