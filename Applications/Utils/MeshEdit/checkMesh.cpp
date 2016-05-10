/**
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <array>
#include <string>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/AABB.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/Legacy/MeshIO.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Checks mesh properties", ' ', BaseLib::BuildInfo::git_describe);
	TCLAP::UnlabeledValueArg<std::string> mesh_arg("mesh-file","input mesh file",true,"","string");
	cmd.add( mesh_arg );
	TCLAP::SwitchArg valid_arg("v","validation","validate the mesh");
	cmd.add( valid_arg );
	TCLAP::SwitchArg print_properties_arg("p","print_properties","print properties stored in the mesh");
	cmd.add( print_properties_arg );

	cmd.parse( argc, argv );

	const std::string filename(mesh_arg.getValue());

	// read the mesh file
	BaseLib::MemWatch mem_watch;
	const unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
	BaseLib::RunTime run_time;
	run_time.start();
	const MeshLib::Mesh* mesh = FileIO::readMeshFromFile(filename); // FileIO outputs nr. of nodes and elements
	if (!mesh)
		return 1;

	const unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
	if (mem_with_mesh>0)
		INFO ("Memory size: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
	INFO ("Time for reading: %g s", run_time.elapsed());

	// Geometric information
	const GeoLib::AABB aabb(MeshLib::MeshInformation::getBoundingBox(*mesh));
	auto minPt(aabb.getMinPoint());
	auto maxPt(aabb.getMaxPoint());
	INFO("Node coordinates:");
	INFO("\tx [%g, %g] (extent %g)", minPt[0], maxPt[0], maxPt[0]-minPt[0]);
	INFO("\ty [%g, %g] (extent %g)", minPt[1], maxPt[1], maxPt[1]-minPt[1]);
	INFO("\tz [%g, %g] (extent %g)", minPt[2], maxPt[2], maxPt[2]-minPt[2]);

	INFO("Edge length: [%g, %g]", mesh->getMinEdgeLength(), mesh->getMaxEdgeLength());

	// Element information
	const std::array<unsigned, 7> nr_ele_types(MeshLib::MeshInformation::getNumberOfElementTypes(*mesh));
	INFO("Element types:");
	unsigned etype = 0;
	if (nr_ele_types[etype]>0) INFO("\t%d lines", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d triangles", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d quads", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d tetrahedra", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d hexahedra", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d pyramids", nr_ele_types[etype]);
	if (nr_ele_types[++etype]>0) INFO("\t%d prisms", nr_ele_types[etype]);

	std::vector<std::string> const& vec_names (mesh->getProperties().getPropertyVectorNames());
	INFO("There are %d properties in the mesh:", vec_names.size());
	for (std::size_t i=0; i<vec_names.size(); ++i)
	{
		auto vec_bounds (MeshLib::MeshInformation::getValueBounds<int>(*mesh, vec_names[i]));
		if (vec_bounds.second != std::numeric_limits<int>::max())
			INFO("\t%s: [%d, %d]", vec_names[i].c_str(), vec_bounds.first, vec_bounds.second)
		else
		{
			auto vec_bounds (MeshLib::MeshInformation::getValueBounds<double>(*mesh, vec_names[i]));
			if (vec_bounds.second != std::numeric_limits<double>::max())
				INFO("\t%s: [%g, %g]", vec_names[i].c_str(), vec_bounds.first, vec_bounds.second)
			else
				INFO("\t%s: Unknown properties", vec_names[i].c_str())
		}
	}

	if (valid_arg.isSet()) {
		// Validation
		MeshLib::MeshValidation validation(*const_cast<MeshLib::Mesh*>(mesh)); // MeshValidation outputs error messages
		/* Remark: MeshValidation can modify the original mesh */

		unsigned const n_holes (MeshLib::MeshValidation::detectHoles(*mesh));
		if (n_holes > 0)
		{
			INFO("%d hole(s) detected within the mesh", n_holes);
		}
		else
			INFO ("No holes found within the mesh.");
	}

	delete mesh;
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();
}
