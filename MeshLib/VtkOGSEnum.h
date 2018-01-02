/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshEnums.h"

MeshLib::CellType VtkCellTypeToOGS(int type);

int OGSToVtkCellType(MeshLib::CellType ogs);
