/**
 * @file ConvertRasterToMesh.h
 * @author Thomas Fischer
 * @date Nov 14, 2012
 *
 * @copyright
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef CONVERTRASTERTOMESH_H_
#define CONVERTRASTERTOMESH_H_

// GeoLib
#include "Raster.h"

// MeshLib
#include "MshEnums.h"

namespace MeshLib {

// forward declaration
class Mesh;

/**
 * Struct gives a selection of possible interpretations for intensities
 */
struct UseIntensityAs
{
	enum type {
		ELEVATION,
		MATERIAL,
		NONE
	};
};

/**
 * Class to convert raster data into meshes. Based on Karsten Rinks algorithm.
 */
class ConvertRasterToMesh {
public:
	ConvertRasterToMesh(GeoLib::Raster const& raster, MshElemType::type elem_type,
					  UseIntensityAs::type intensity_type);
	~ConvertRasterToMesh();
	MeshLib::Mesh* execute() const;
private:
	double getExistingValue(GeoLib::Raster::const_iterator beg, GeoLib::Raster::const_iterator last) const;
	MeshLib::Mesh* constructMesh(const double* pix_vals, const bool* vis_nodes) const;
	GeoLib::Raster const& _raster;
	MshElemType::type _elem_type;
	UseIntensityAs::type _intensity_type;
};

} // end namespace MeshLib

#endif /* CONVERTRASTERTOMESH_H_ */
