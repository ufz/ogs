/**
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "FileIO/readMeshFromFile.h"
#include "FileIO/Legacy/MeshIO.h"
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Checks mesh properties", ' ', BaseLib::BuildInfo::git_version_sha1);
	TCLAP::ValueArg<std::string> mesh_arg("i","input","input mesh file",true,"","string");
	cmd.add( mesh_arg );

	cmd.parse( argc, argv );

	const std::string filename(mesh_arg.getValue());

	// read the mesh file
	BaseLib::MemWatch mem_watch;
	const unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
	BaseLib::RunTime run_time;
	run_time.start();
	const MeshLib::Mesh* mesh = FileIO::readMeshFromFile(filename); // FileIO outputs nr. of nodes and elements
	run_time.stop();
	if (!mesh)
		return 1;

	const unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
	if (mem_with_mesh>0)
		INFO ("Memory size: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
	INFO ("Time for reading: %g s", run_time.elapsed());

	// Geometric information
	const GeoLib::AABB<MeshLib::Node> aabb(MeshLib::MeshInformation::getBoundingBox(*mesh));
	auto minPt(aabb.getMinPoint());
	auto maxPt(aabb.getMaxPoint());
	INFO("Node coordinates:");
	INFO("\tx [%g, %g] range %g", minPt[0], maxPt[0], maxPt[0]-minPt[0]);
	INFO("\ty [%g, %g] range %g", minPt[1], maxPt[1], maxPt[1]-minPt[1]);
	INFO("\tz [%g, %g] range %g", minPt[2], maxPt[2], maxPt[2]-minPt[2]);

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

	const std::pair<unsigned, unsigned> minmax_values(MeshLib::MeshInformation::getValueBounds(*mesh));
	INFO("Material IDs: [%d, %d]", minmax_values.first, minmax_values.second);


	delete mesh;
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();
}
