/**
 * \file DirectConditionGenerator.cpp
 * 2012/01/04 KR Initial implementation
 *
 */

#include "DirectConditionGenerator.h"

#include "OGSRaster.h"
#include "MshEditor.h"
#include "PointWithID.h"


const std::vector< std::pair<size_t,double> > DirectConditionGenerator::fromRasterToSurfaceNodes(const MeshLib::CFEMesh &mesh, const std::string &filename)
{
	if (_direct_values.empty())
	{
		double origin_x(0), origin_y(0), delta(0);
		size_t imgwidth(0), imgheight(0);

		double* img = OGSRaster::loadDataFromASC(filename, origin_x, origin_y, imgwidth, imgheight, delta);

		const std::vector<GEOLIB::PointWithID*> surface_nodes ( MshEditor::getSurfaceNodes(mesh) );		
		//std::vector<MeshLib::CNode*> nodes = mesh.nod_vector;
		const size_t nNodes(surface_nodes.size());
		for (size_t i=0; i<nNodes; i++)
		{
			const double* coords (surface_nodes[i]->getData());

			if (coords[0]>=origin_x && coords[0]<(origin_x+(delta*imgwidth)) && coords[1]>=origin_y && coords[1]<(origin_y+(delta*imgheight)))
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

