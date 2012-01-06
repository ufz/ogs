/**
 * \file DirectConditionGenerator.cpp
 * 2012/01/04 KR Initial implementation
 *
 */

#include "DirectConditionGenerator.h"

#include "OGSRaster.h"
#include "PointWithID.h"


const std::vector< std::pair<size_t,double> > DirectConditionGenerator::fromRasterToSurfaceNodes(const MeshLib::CFEMesh &mesh, const std::string &filename)
{
	if (_direct_values.empty())
	{
		double origin_x(0), origin_y(0), delta(0);
		size_t imgwidth(0), imgheight(0);

		double* img = OGSRaster::loadDataFromASC(filename, origin_x, origin_y, imgwidth, imgheight, delta);

		const std::vector<GEOLIB::PointWithID*> surface_nodes ( this->getSurfaceNodes(mesh) );		
		//std::vector<MeshLib::CNode*> nodes = mesh.nod_vector;
		const size_t nNodes(surface_nodes.size());
		for (size_t i=0; i<nNodes; i++)
		{
			const double* coords (surface_nodes[i]->getData());

			if (coords[0]>=origin_x && coords[0]<(origin_x+(delta*imgheight)) && coords[1]>=origin_y && coords[1]<(origin_y+(delta*imgwidth)))
			{
				size_t cell_x = static_cast<size_t>(floor((coords[0] - origin_x)/delta));
				size_t cell_y = static_cast<size_t>(floor((coords[1] - origin_y)/delta));
				size_t index = cell_y*imgwidth+cell_x;
				if (img[index] != -9999)
					_direct_values.push_back( std::pair<size_t, double>(surface_nodes[i]->getID(),img[index]) );
			}
		}

		delete[] img;

	}
	else
		std::cout << "Error in DiretConditionGenerator::fromRasterToNodes() - Data vector contains outdated values..." << std::endl;

	return _direct_values;
}


int DirectConditionGenerator::writeToFile(const std::string &name) const
{
	std::ofstream out( name.c_str(), std::ios::out );
	
	if (out)
	{
		for (std::vector< std::pair<size_t,double> >::const_iterator it = _direct_values.begin(); it != _direct_values.end(); ++it)
			out << it->first << "\t" << it->second << std::endl;

		out.close();
	}
	return 0;
}

const std::vector<GEOLIB::PointWithID*> DirectConditionGenerator::getSurfaceNodes(const MeshLib::CFEMesh &mesh)
{
	// Sort points lexicographically
	size_t nNodes (mesh.nod_vector.size());
	std::vector<GEOLIB::PointWithID*> nodes;
	std::vector<size_t> perm;
	for (size_t j(0); j<nNodes; j++)
	{
		nodes.push_back(new GEOLIB::PointWithID(mesh.nod_vector[j]->getData(), j));		
		perm.push_back(j);
	}
	Quicksort<GEOLIB::PointWithID*> (nodes, 0, nodes.size(), perm);

	// Extract surface points
	double eps (sqrt(std::numeric_limits<double>::min()));
	std::vector<GEOLIB::PointWithID*> surface_pnts;
	for (size_t k(1); k < nNodes; k++)
	{
		const GEOLIB::PointWithID& p0 (*(nodes[k - 1]));
		const GEOLIB::PointWithID& p1 (*(nodes[k]));
		if (fabs (p0[0] - p1[0]) > eps || fabs (p0[1] - p1[1]) > eps)
			surface_pnts.push_back (nodes[k - 1]);
	}
	// Add last point
	surface_pnts.push_back (nodes[nNodes - 1]);
	return surface_pnts;
}