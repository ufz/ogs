/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-04
 * \brief  Implementation of the DirectConditionGenerator class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <fstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "DirectConditionGenerator.h"

#include "Raster.h"
#include "MeshSurfaceExtraction.h"
#include "PointWithID.h"
#include "Mesh.h"
#include "Node.h"

#include <cmath>
#include <limits>

const std::vector< std::pair<size_t,double> >& DirectConditionGenerator::directToSurfaceNodes(const MeshLib::Mesh &mesh, const std::string &filename)
{
	if (_direct_values.empty())
	{
		GeoLib::Raster* raster (GeoLib::Raster::readRaster(filename));
		if (! raster) {
			ERR("Error in DirectConditionGenerator::directToSurfaceNodes() - could not load raster file.");
			return _direct_values;
		}

		const MathLib::Vector3 dir(0,0,-1);
		const std::vector<GeoLib::PointWithID*> surface_nodes(MeshLib::MeshSurfaceExtraction::getSurfaceNodes(mesh, dir) );
		const size_t nNodes(surface_nodes.size());
		const double no_data (raster->getNoDataValue());
		_direct_values.reserve(nNodes);
		for (size_t i=0; i<nNodes; i++)
		{
			double val (raster->getValueAtPoint(*surface_nodes[i]));
			val = (fabs(val-no_data) < std::numeric_limits<double>::epsilon()) ? 0 : val;
			_direct_values.push_back (std::pair<size_t, double>(surface_nodes[i]->getID(), val));
		}
		delete raster;
	}
	else
		ERR("Error in DirectConditionGenerator::directToSurfaceNodes() - Data vector contains outdated values.");

	return _direct_values;
}


const std::vector< std::pair<size_t,double> >& DirectConditionGenerator::directWithSurfaceIntegration(MeshLib::Mesh &mesh, const std::string &filename, double scaling)
{
	if (_direct_values.empty())
	{
		GeoLib::Raster* raster (GeoLib::Raster::readRaster(filename));
		if (!raster) {
			ERR("Error in DirectConditionGenerator::directWithSurfaceIntegration() - could not load raster file.");
			return _direct_values;
		}

		const MathLib::Vector3 dir(0,0,-1);
		MeshLib::Mesh* sfc_mesh (MeshLib::MeshSurfaceExtraction::getMeshSurface(mesh, dir, true));
		std::vector<double> node_area_vec;
		MeshLib::MeshSurfaceExtraction::getSurfaceAreaForNodes(*sfc_mesh, node_area_vec);
		const std::vector<MeshLib::Node*> &surface_nodes (sfc_mesh->getNodes());
		const size_t nNodes(sfc_mesh->getNNodes());
		const double no_data (raster->getNoDataValue());
		_direct_values.reserve(nNodes);
		for (size_t i=0; i<nNodes; ++i)
		{
			double val (raster->getValueAtPoint(*surface_nodes[i]));
			val = (fabs(val-no_data) < std::numeric_limits<double>::epsilon()) ? 0 : (val*node_area_vec[i]*scaling);
			_direct_values.push_back (std::pair<size_t, double>(surface_nodes[i]->getID(), val));
		}

		delete raster;
	}
	else
		std::cout << "Error in DirectConditionGenerator::directWithSurfaceIntegration() - Data vector contains outdated values..." << std::endl;

	return _direct_values;
}


int DirectConditionGenerator::writeToFile(const std::string &name) const
{
	std::ofstream out( name.c_str(), std::ios::out );

	if (out)
	{
		for (std::vector< std::pair<size_t,double> >::const_iterator it = _direct_values.begin(); it != _direct_values.end(); ++it)
			out << it->first << "\t" << it->second << "\n";

		out.close();
	}
	return 0;
}

