/**
 * @file ConvertRasterToMesh.cpp
 * @author Thomas Fischer
 * @date Nov 14, 2012
 * @brief Implementation of the ConvertRasterToMesh class.
 *
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConvertRasterToMesh.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Elements/Quad.h"

namespace MeshLib {

ConvertRasterToMesh::ConvertRasterToMesh(GeoLib::Raster const& raster, MeshElemType elem_type,
				UseIntensityAs intensity_type) :
	_raster(raster), _elem_type(elem_type), _intensity_type(intensity_type)
{}

ConvertRasterToMesh::~ConvertRasterToMesh()
{}

MeshLib::Mesh* ConvertRasterToMesh::execute() const
{
	GeoLib::RasterHeader const& header (_raster.getHeader());
	const std::size_t height(header.n_rows+1);
	const std::size_t width(header.n_cols+1);
	double* pix_vals(new double[height*width]);
	bool* vis_nodes(new bool[height*width]);

	// determine a valid value for substitution of no data values
	double substitution(getExistingValue(_raster.begin(), _raster.end()));

	// fill first row with non visual nodes
	for (std::size_t j = 0; j < width; j++) {
		pix_vals[j] = 0;
		vis_nodes[j] = false;
	}

	GeoLib::Raster::const_iterator raster_it(_raster.begin());
	for (std::size_t i = 0; i < height; ++i) {
		for (std::size_t j = 0; j < width; ++j) {
			const std::size_t index = (i+1) * width + j;
			if (*raster_it == header.no_data) {
				pix_vals[index] = substitution;
				vis_nodes[index] = false;
			} else {
				pix_vals[index] = *raster_it;
				vis_nodes[index] = true;
			}
			++raster_it;
		}
		// fill last column with non-visual nodes
		pix_vals[(i + 2) * width - 1] = 0;
		vis_nodes[(i + 2) * width - 1] = false;
	}

	MeshLib::Mesh* mesh = constructMesh(pix_vals, vis_nodes);

	delete [] pix_vals;
	delete [] vis_nodes;

	return mesh;
}

MeshLib::Mesh* ConvertRasterToMesh::constructMesh(const double* pix_vals, const bool* vis_nodes) const
{
	GeoLib::RasterHeader const& header (_raster.getHeader());
	const std::size_t height (header.n_rows+1);
	const std::size_t width (header.n_cols+1);
	std::size_t node_idx_count(0);
	const double distance(header.cell_size);

	const std::size_t size(height*width);
	int* node_idx_map(new int[size]);
	for (std::size_t k(0); k<size; ++k) node_idx_map[k] = -1;

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	for (std::size_t i = 0; i < height; i++) {
		for (std::size_t j = 0; j < width; j++) {
			const std::size_t index = i * width + j;
			double zValue = (_intensity_type == UseIntensityAs::ELEVATION) ? 
				pix_vals[index] : header.origin[2];
			MeshLib::Node* node(new MeshLib::Node(header.origin[0] + (distance * j), 
			                                      header.origin[1] + (distance * i), zValue));
			nodes.push_back(node);
			node_idx_map[index] = node_idx_count;
			node_idx_count++;
		}
	}

	std::vector<int> mat_ids;
	// set mesh elements
	for (std::size_t i = 0; i < height; i++) {
		for (std::size_t j = 0; j < width; j++) {
			const int index = i * width + j;
			if ((node_idx_map[index] != -1) && (node_idx_map[index + 1] != -1)
					&& (node_idx_map[index + width] != -1)
					&& (node_idx_map[index + width + 1] != -1)
					&& (vis_nodes[index + width])) {
				const int mat = (_intensity_type != UseIntensityAs::DATAVECTOR) ? 0
								: static_cast<int> (pix_vals[index + width]);
				if (_elem_type == MeshElemType::TRIANGLE) {
					MeshLib::Node** tri1_nodes = new MeshLib::Node*[3];
					tri1_nodes[0] = nodes[node_idx_map[index]];
					tri1_nodes[1] = nodes[node_idx_map[index + 1]];
					tri1_nodes[2] = nodes[node_idx_map[index + width]];

					MeshLib::Node** tri2_nodes = new MeshLib::Node*[3];
					tri2_nodes[0] = nodes[node_idx_map[index + 1]];
					tri2_nodes[1] = nodes[node_idx_map[index + width + 1]];
					tri2_nodes[2] = nodes[node_idx_map[index + width]];

					// upper left triangle
					elements.push_back(new MeshLib::Tri(tri1_nodes));
					mat_ids.push_back(mat);
					// lower right triangle
					elements.push_back(new MeshLib::Tri(tri2_nodes));
					mat_ids.push_back(mat);
				}
				if (_elem_type == MeshElemType::QUAD) {
					MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
					quad_nodes[0] = nodes[node_idx_map[index]];
					quad_nodes[1] = nodes[node_idx_map[index + 1]];
					quad_nodes[2] = nodes[node_idx_map[index + width + 1]];
					quad_nodes[3] = nodes[node_idx_map[index + width]];
					elements.push_back(new MeshLib::Quad(quad_nodes));
					mat_ids.push_back(mat);
				}
			}
		}
	}
	delete [] node_idx_map;

	// the name is only a temp-name, the name given in the dialog is set later
	MeshLib::Properties properties;
	boost::optional<MeshLib::PropertyVector<int>&> materials =
	    properties.createNewPropertyVector<int>("MaterialIDs",
	                                            MeshLib::MeshItemType::Cell);
	assert(materials != boost::none);
	materials->reserve(mat_ids.size());
	std::copy(mat_ids.cbegin(), mat_ids.cend(), std::back_inserter(*materials));
	return new MeshLib::Mesh("RasterDataMesh", nodes, elements, properties);
}


double ConvertRasterToMesh::getExistingValue(GeoLib::Raster::const_iterator first, GeoLib::Raster::const_iterator last) const
{
	double const no_data (_raster.getHeader().no_data);
	for (GeoLib::Raster::const_iterator it(first); it != last; ++it) {
		if (*it != no_data)
			return *it;
	}
	return no_data;
}
} // end namespace MeshLib
