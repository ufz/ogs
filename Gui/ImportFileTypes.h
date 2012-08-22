/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file ImportFileTypes.h
 *
 * Created on 2012-08-20 by Karsten Rink
 */

#ifndef IMPORTFILETYPES_H
#define IMPORTFILETYPES_H

#include <string>

/**
 * \brief Types of supported import file formats.
 */
struct ImportFileType
{
	enum type {
		INVALID = 0,
		GMS,
		GMSH,
		NETCDF,
		OGS,
		PETREL,
		RASTER,
		SHAPE,
		TETGEN,
		VTK
	};
};

#endif //IMPORTFILETYPES_H
