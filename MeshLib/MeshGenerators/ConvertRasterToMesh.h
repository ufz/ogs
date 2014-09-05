/**
 * @file ConvertRasterToMesh.h
 * @author Thomas Fischer
 * @date Nov 14, 2012
 * @brief Definition of the ConvertRasterToMesh class.
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CONVERTRASTERTOMESH_H_
#define CONVERTRASTERTOMESH_H_

// GeoLib
#include "Raster.h"

// MeshLib
#include "MeshEnums.h"

namespace MeshLib {

// forward declaration
class Mesh;

/**
 * Struct gives a selection of possible interpretations for intensities
 */
enum class UseIntensityAs
{
	ELEVATION,
	MATERIAL,
	NONE
};

/**
 * Class to convert raster data into meshes. Based on Karsten Rinks algorithm.
 */
class ConvertRasterToMesh {
public:
	ConvertRasterToMesh(GeoLib::Raster const& raster, MeshElemType elem_type,
					  UseIntensityAs intensity_type);
	~ConvertRasterToMesh();
	MeshLib::Mesh* execute() const;
private:
	double getExistingValue(GeoLib::Raster::const_iterator beg, GeoLib::Raster::const_iterator last) const;
	MeshLib::Mesh* constructMesh(const double* pix_vals, const bool* vis_nodes) const;
	GeoLib::Raster const& _raster;
	MeshElemType _elem_type;
	UseIntensityAs _intensity_type;
};

} // end namespace MeshLib

#endif /* CONVERTRASTERTOMESH_H_ */
