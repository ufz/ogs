/**
 * @file ConvertRasterToMesh.cpp
 * @author Thomas Fischer
 * @date Nov 14, 2012
 * @brief Implementation of the ConvertRasterToMesh class.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "ConvertRasterToMesh.h"

// MeshLib
#include "Node.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"

namespace MeshLib {

ConvertRasterToMesh::ConvertRasterToMesh(GeoLib::Raster const& raster, MeshElemType elem_type,
				UseIntensityAs intensity_type) :
	_raster(raster), _elem_type(elem_type), _intensity_type(intensity_type)
{}

ConvertRasterToMesh::~ConvertRasterToMesh()
{}

MeshLib::Mesh* ConvertRasterToMesh::execute() const
{
	const size_t height(_raster.getNRows()+1);
	const size_t width(_raster.getNCols()+1);
	const size_t size(height*width);
	double* pix_vals(new double[size]);
	bool* vis_nodes(new bool[size]);

	// determine a valid value for substitution of no data values
	double substitution(getExistingValue(_raster.begin(), _raster.end()));

	// fill first row with non visual nodes
	for (size_t j = 0; j < _raster.getNCols(); j++) {
		pix_vals[j] = 0;
		vis_nodes[j] = false;
	}

	GeoLib::Raster::const_iterator raster_it(_raster.begin());
	for (size_t i = 0; i < _raster.getNRows(); ++i) {
		for (size_t j = 0; j < _raster.getNCols(); ++j) {
			const size_t index = (i+1) * width + j;
			if (*raster_it == _raster.getNoDataValue()) {
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
	const size_t height = _raster.getNRows()+1;
	const size_t width = _raster.getNCols()+1;
	size_t node_idx_count(0);
	const double distance(_raster.getRasterPixelDistance());
	const double x_offset(_raster.getOrigin()[0]); // - distance / 2.0);
	const double y_offset(_raster.getOrigin()[1]); // - distance / 2.0);

	const size_t size(height*width);
	int* node_idx_map(new int[size]);
	for (std::size_t k(0); k<size; ++k) node_idx_map[k] = -1;

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	for (size_t i = 0; i < height; i++) {
		for (size_t j = 0; j < width; j++) {
			const size_t index = i * width + j;

//			bool set_node(true);
//			bool set_node(false);
//			if (j == 0 && i == height)
//				set_node = vis_nodes[index];
//			else if (j == 0)
//				set_node = (vis_nodes[index] || vis_nodes[index + height]);
//			else if (i == width)
//				set_node = (vis_nodes[index] || vis_nodes[index - 1]);
//			else set_node = (vis_nodes[index] || vis_nodes[index - 1] || vis_nodes[index + height]
//							|| vis_nodes[index + height - 1]);
//			if (set_node) {
				double zValue = (_intensity_type == UseIntensityAs::ELEVATION) ? pix_vals[index]
								: _raster.getOrigin()[2];
				MeshLib::Node* node(new MeshLib::Node(x_offset + (distance * j), y_offset
								+ (distance * i), zValue));
				nodes.push_back(node);
				node_idx_map[index] = node_idx_count;
				node_idx_count++;
//			}
		}
	}

	// set mesh elements
	for (size_t i = 0; i < _raster.getNRows(); i++) {
		for (size_t j = 0; j < _raster.getNCols(); j++) {
			const int index = i * width + j;
			if ((node_idx_map[index] != -1) && (node_idx_map[index + 1] != -1)
					&& (node_idx_map[index + width] != -1)
					&& (node_idx_map[index + width + 1] != -1)
					&& (vis_nodes[index + width])) {
				const int mat = (_intensity_type != UseIntensityAs::MATERIAL) ? 0
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

					elements.push_back(new MeshLib::Tri(tri1_nodes, mat)); // upper left triangle
					elements.push_back(new MeshLib::Tri(tri2_nodes, mat)); // lower right triangle
				}
				if (_elem_type == MeshElemType::QUAD) {
					MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
					quad_nodes[0] = nodes[node_idx_map[index]];
					quad_nodes[1] = nodes[node_idx_map[index + 1]];
					quad_nodes[2] = nodes[node_idx_map[index + width + 1]];
					quad_nodes[3] = nodes[node_idx_map[index + width]];
					elements.push_back(new MeshLib::Quad(quad_nodes, mat));
				}
			}
		}
	}
	delete [] node_idx_map;

	return new MeshLib::Mesh("RasterDataMesh", nodes, elements); // the name is only a temp-name, the name given in the dialog is set later
}


double ConvertRasterToMesh::getExistingValue(GeoLib::Raster::const_iterator first, GeoLib::Raster::const_iterator last) const
{
	for (GeoLib::Raster::const_iterator it(first); it != last; ++it) {
		if (*it != _raster.getNoDataValue())
			return *it;
	}
	return _raster.getNoDataValue();
}
} // end namespace MeshLib
