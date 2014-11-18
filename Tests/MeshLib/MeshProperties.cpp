/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *			Distributed under a Modified BSD License.
 *			  See accompanying file LICENSE.txt or
 *			  http://www.opengeosys.org/LICENSE.txt
 */

#include <ctime>
#include "gtest/gtest.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Line.h"

class MeshLibMeshProperties : public ::testing::Test
{
	public:
	MeshLibMeshProperties()
		: mesh(nullptr)
	{
		mesh = MeshLib::MeshGenerator::generateRegularHexMesh(1.0, mesh_size);
	}

	~MeshLibMeshProperties()
	{
		delete mesh;
	}

	static std::size_t const mesh_size = 9;
	MeshLib::Mesh * mesh;
};
std::size_t const MeshLibMeshProperties::mesh_size;

TEST_F(MeshLibMeshProperties, AddDoubleProperties)
{
	ASSERT_TRUE(mesh != nullptr);
	const std::size_t size(mesh_size*mesh_size*mesh_size);
	std::vector<double> double_properties(size);
	std::iota(double_properties.begin(), double_properties.end(), 1);

	std::string const& prop_name("FirstTestProperty");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name, double_properties);

	boost::optional<std::vector<double> const&>
		double_properties_cpy(mesh->getProperties().getProperty<double>(prop_name));

	for (std::size_t k(0); k<size; k++) {
		ASSERT_EQ(double_properties[k], (*double_properties_cpy)[k]);
	}
}

