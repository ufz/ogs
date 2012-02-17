/**
 * \file DirectConditionGenerator.cpp
 * 2012/01/04 KR Initial implementation
 *
 */

#include "DirectConditionGenerator.h"

#include "OGSRaster.h"
#include "MshEditor.h"
#include "PointWithID.h"

#include "fem_ele.h"

const std::vector< std::pair<size_t,double> >& DirectConditionGenerator::directToSurfaceNodes(const MeshLib::CFEMesh &mesh, const std::string &filename)
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
		std::cout << "Error in DiretConditionGenerator::directToSurfaceNodes() - Data vector contains outdated values..." << std::endl;

	return _direct_values;
}


const std::vector< std::pair<size_t,double> >& DirectConditionGenerator::directWithSurfaceIntegration(MeshLib::CFEMesh &mesh, const std::string &filename)
{
	double no_data_value = -9999; // TODO: get this from asc-reader!
	double node_val[8] = {0,0,0,0,0,0,0,0}; // maximum possible number of nodes per face (just in case ...)

	if (_direct_values.empty())
	{
		mesh.MarkInterface_mHM_Hydro_3D(); // mark element faces on the surface

		double origin_x(0), origin_y(0), delta(0);
		size_t imgwidth(0), imgheight(0);

		double* img = OGSRaster::loadDataFromASC(filename, origin_x, origin_y, imgwidth, imgheight, delta);

		// copied from CFEMesh::Precipitation2NeumannBC() by WW
		size_t nFaces = mesh.face_vector.size();
		for(size_t i=0; i<nFaces; i++)
		{
			MeshLib::CElem* elem = mesh.face_vector[i];
			if(!elem->GetMark())
				continue;

			size_t nElemNodes = elem->GetNodesNumber(false);
			for(size_t k=0; k<nElemNodes; k++)
				node_val[k] = 0.0;

			for(size_t k=0; k<nElemNodes; k++)
			{
				//MeshLib::CNode* node = elem->GetNode(k);
				double const* const pnt_k (elem->GetNode(k)->getData());

				size_t nx = static_cast<size_t>(floor((pnt_k[0] - origin_x) / delta));
				size_t ny = static_cast<size_t>(floor((pnt_k[1] - origin_y) / delta));
				/*
				if(ny < 0)
					ny = 0;
				if(ny > static_cast<long>(nrows))
					ny = nrows;

				if(nx * csize + x0 >= pnt_k[0])
					nx -= 1;
				if(ny * csize + y0 >= pnt_k[1])
					ny -= 1;
				if(nx >= static_cast<long>(ncols) - 1)
					nx = ncols - 2;
				if(ny >= static_cast<long>(nrows) - 1)
					ny = nrows - 2;
				if(nx < 0)
					nx = 0;
				if(ny < 0)
					ny = 0;

				node_val[k] = zz[ncols * ny + nx];
				*/
				node_val[k] = img[ny * imgwidth + nx];
				if (fabs(node_val[k] - no_data_value) < std::numeric_limits<double>::min())
					node_val[k] = 0.;
			}

			elem->ComputeVolume();

			FiniteElement::CElement* fem ( NULL );
			fem->setOrder(mesh.getOrder() + 1);
			fem->ConfigElement(elem);
			fem->FaceIntegration(node_val);

			for(size_t k=0; k<elem->GetNodesNumber(false); k++)
			{
				MeshLib::CNode* node = elem->GetNode(k);
				node->SetMark(true);
				//val[node->GetIndex()] += node_val[k];
				_direct_values.push_back( std::pair<size_t, double>(node->GetIndex(), node_val[k]) );
			}
		}


	}
	else
		std::cout << "Error in DiretConditionGenerator::directWithSurfaceIntegration() - Data vector contains outdated values..." << std::endl;

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

