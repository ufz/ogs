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
#include "Location.h"
#include "PropertyVector.h"

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

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

	template <typename T>
	MeshLib::PropertyVector<T*>* createGroupPropertyValueVector(
		std::size_t n_mat_groups)
	{
		const std::size_t n_elements(mesh_size*mesh_size*mesh_size);
		std::vector<std::size_t> mat_group_idx_map(n_elements);
		// create simple mat_group to index mapping
		for (std::size_t j(0); j<n_mat_groups; j++) {
			std::size_t const lower((double)(j)/(double)(n_mat_groups)*n_elements);
			std::size_t const upper((double)(j+1)/(double)(n_mat_groups)*n_elements);
			for (std::size_t k(lower); k<upper; k++) {
				mat_group_idx_map[k] = j;
			}
		}
		MeshLib::PropertyVector<T*> *group_props(
			new MeshLib::PropertyVector<T*>(
				n_mat_groups, mat_group_idx_map
			)
		);
		for (auto it=group_props->begin(); it != group_props->end(); it++) {
			(*it) = new T;
			for (std::size_t idx(0); idx<(*it)->size(); idx++) {
				(*(*it))[idx] = std::distance(group_props->begin(), it)+idx;
			}
		}

		return group_props;
	}

	static std::size_t const mesh_size = 5;
	MeshLib::Mesh * mesh;
};
std::size_t const MeshLibMeshProperties::mesh_size;

TEST_F(MeshLibMeshProperties, AddDoubleProperties)
{
	ASSERT_TRUE(mesh != nullptr);
	const std::size_t size(mesh_size*mesh_size*mesh_size);
	MeshLib::PropertyVector<double> *double_properties(
		new MeshLib::PropertyVector<double>(size)
	);
	// init property values
	std::iota(double_properties->begin(), double_properties->end(), 1);

	std::string const& prop_name("FirstTestProperty");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name, double_properties,
		MeshLib::MeshItemType::Cell);

	boost::optional<MeshLib::PropertyVector<double> *>
		double_properties_cpy(mesh->getProperties().getProperty<double>(
			prop_name, MeshLib::MeshItemType::Cell
		));
	ASSERT_FALSE(!double_properties_cpy);

	for (std::size_t k(0); k<size; k++) {
		ASSERT_EQ((*double_properties)[k], (*(*double_properties_cpy))[k]);
	}

	mesh->getProperties().removeProperty(prop_name, MeshLib::MeshItemType::Cell);
	boost::optional<MeshLib::PropertyVector<double> *>
		removed_double_properties(mesh->getProperties().getProperty<double>(prop_name,
			MeshLib::MeshItemType::Cell)
		);

	ASSERT_TRUE(!removed_double_properties);
}

TEST_F(MeshLibMeshProperties, AddDoublePointerProperties)
{
	ASSERT_TRUE(mesh != nullptr);
	const std::size_t n_mat_groups(10);
	const std::size_t n_elements(mesh_size*mesh_size*mesh_size);
	std::vector<std::size_t> mat_group_idx_map(n_elements);
	// create simple mat_group to index mapping
	for (std::size_t j(0); j<n_mat_groups; j++) {
		std::size_t const lower((double)(j)/(double)(n_mat_groups)*n_elements);
		std::size_t const upper((double)(j+1)/(double)(n_mat_groups)*n_elements);
		for (std::size_t k(lower); k<upper; k++) {
			mat_group_idx_map[k] = j;
		}
	}
	MeshLib::PropertyVector<double*> *pointer_properties(
		new MeshLib::PropertyVector<double*>(n_mat_groups, mat_group_idx_map)
	);
	for (auto it=pointer_properties->begin(); it != pointer_properties->end(); it++)
	{
		(*it) = new double;
		*(*it) = std::distance(pointer_properties->begin(), it);
	}

	std::string const& prop_name("TestMatGroupProperty");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name, pointer_properties,
		MeshLib::MeshItemType::Cell);

	boost::optional<MeshLib::PropertyVector<double*> *>
		pointer_properties_cpy(mesh->getProperties().getProperty<double*>(
			prop_name, MeshLib::MeshItemType::Cell
		));
	ASSERT_FALSE(!pointer_properties_cpy);

	for (std::size_t k(0); k<n_elements; k++) {
		ASSERT_EQ((*pointer_properties)[k], (*(*pointer_properties_cpy))[k]);
	}

	mesh->getProperties().removeProperty(prop_name, MeshLib::MeshItemType::Cell);
	boost::optional<MeshLib::PropertyVector<double*> *>
		removed_double_properties(mesh->getProperties().getProperty<double*>(
			prop_name, MeshLib::MeshItemType::Cell
		));

	ASSERT_TRUE(!removed_double_properties);
}

TEST_F(MeshLibMeshProperties, AddArrayPointerProperties)
{
	ASSERT_TRUE(mesh != nullptr);
	const std::size_t n_mat_groups(10);
	const std::size_t n_elements(mesh_size*mesh_size*mesh_size);
	std::vector<std::size_t> mat_group_idx_map(n_elements);
	// create simple mat_group to index mapping
	for (std::size_t j(0); j<n_mat_groups; j++) {
		std::size_t const lower((double)(j)/(double)(n_mat_groups)*n_elements);
		std::size_t const upper((double)(j+1)/(double)(n_mat_groups)*n_elements);
		for (std::size_t k(lower); k<upper; k++) {
			mat_group_idx_map[k] = j;
		}
	}
	MeshLib::PropertyVector<std::array<double,3>*> *pointer_properties(
		new MeshLib::PropertyVector<std::array<double,3>*>(
			n_mat_groups, mat_group_idx_map
		)
	);
	for (auto it=pointer_properties->begin(); it != pointer_properties->end(); it++)
	{
		(*it) = new std::array<double,3>;
		(*(*it))[0] = std::distance(pointer_properties->begin(), it);
		(*(*it))[1] = std::distance(pointer_properties->begin(), it)+1;
		(*(*it))[2] = std::distance(pointer_properties->begin(), it)+2;
	}

	std::string const& prop_name("TestMatGroupProperty");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name, pointer_properties,
		MeshLib::MeshItemType::Cell);

	boost::optional<MeshLib::PropertyVector<std::array<double,3>*> *>
		pointer_properties_cpy(
			mesh->getProperties().getProperty<std::array<double,3>*>(
				prop_name, MeshLib::MeshItemType::Cell
			)
		);
	ASSERT_FALSE(!pointer_properties_cpy);

	for (std::size_t k(0); k<n_elements; k++) {
		ASSERT_EQ((*((*pointer_properties)[k]))[0],
			(*((*(*pointer_properties_cpy))[k]))[0]);
		ASSERT_EQ((*((*pointer_properties)[k]))[1],
			(*((*(*pointer_properties_cpy))[k]))[1]);
		ASSERT_EQ((*((*pointer_properties)[k]))[2],
			(*((*(*pointer_properties_cpy))[k]))[2]);
	}

	mesh->getProperties().removeProperty(prop_name, MeshLib::MeshItemType::Cell);
	boost::optional<MeshLib::PropertyVector<std::array<double, 3>*> *>
		removed_double_properties(
			mesh->getProperties().getProperty<std::array<double,3>*>(
				prop_name, MeshLib::MeshItemType::Cell
			)
		);

	ASSERT_TRUE(!removed_double_properties);
}

TEST_F(MeshLibMeshProperties, AddVariousDifferentProperties)
{
	ASSERT_TRUE(mesh != nullptr);

	const std::size_t n_mat_groups(10);
	std::string const& prop_name("TestMatGroupVectorProperty");
	// add a vector with property values to the mesh
	MeshLib::PropertyVector<std::array<double,3>*> *pointer_properties(
		createGroupPropertyValueVector<std::array<double,3>>(n_mat_groups));
	mesh->getProperties().addProperty(prop_name, pointer_properties,
		MeshLib::MeshItemType::Cell);

	// fetch the vector filled with property values from mesh
	boost::optional<MeshLib::PropertyVector<std::array<double,3>*> *>
		pointer_properties_cpy(
			mesh->getProperties().getProperty<std::array<double,3>*>(
				prop_name, MeshLib::MeshItemType::Cell
			)
		);
	ASSERT_FALSE(!pointer_properties_cpy);

	// compare the content
	const std::size_t n_elements(mesh_size*mesh_size*mesh_size);
	for (std::size_t k(0); k<n_elements; k++) {
		ASSERT_EQ((*((*pointer_properties)[k]))[0],
			(*((*(*pointer_properties_cpy))[k]))[0]);
		ASSERT_EQ((*((*pointer_properties)[k]))[1],
			(*((*(*pointer_properties_cpy))[k]))[1]);
		ASSERT_EQ((*((*pointer_properties)[k]))[2],
			(*((*(*pointer_properties_cpy))[k]))[2]);
	}

	// add a 2nd property
	const std::size_t size(mesh_size*mesh_size*mesh_size);
	MeshLib::PropertyVector<std::array<float,9>> *matrix_properties(
		new MeshLib::PropertyVector<std::array<float,9>>(size)
	);
	// init property values
	for (auto it=matrix_properties->begin(); it != matrix_properties->end(); it++)
	{
		for (std::size_t k(0); k<it->size(); k++) {
			(*it)[k] = std::distance(matrix_properties->begin(), it)+k;
		}
	}

	std::string const& prop_name_2("MatrixProperties");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name_2, matrix_properties,
		MeshLib::MeshItemType::Cell);

	// fetch the vector in order to compare the content
	boost::optional<MeshLib::PropertyVector<std::array<float,9>> *>
		matrix_properties_cpy(mesh->getProperties().getProperty<std::array<float,9>>(
			prop_name_2, MeshLib::MeshItemType::Cell
		)
	);
	ASSERT_FALSE(!matrix_properties_cpy);

	// compare the values/matrices
	for (std::size_t k(0); k<size; k++) {
		ASSERT_EQ((*matrix_properties)[k], (*(*matrix_properties_cpy))[k]);
	}

	// add a 3rd property
#ifdef OGS_USE_EIGEN
	MeshLib::PropertyVector<Eigen::Matrix<double,3,3,Eigen::RowMajor>>
		*properties_3(
			new MeshLib::PropertyVector<Eigen::Matrix<double,3,3,Eigen::RowMajor>>(size)
	);

	// init property values
	for (auto it=properties_3->begin(); it != properties_3->end(); it++)
	{
		for (int r(0); r<it->rows(); r++) {
			for (int c(0); c<it->cols(); c++) {
				(*it)(r,c) = std::distance(properties_3->begin(),it)
					+ r*it->cols()+c+1;
			}
		}
	}

	std::string const& prop_name_3("Properties3");
	// add a vector with property values to the mesh
	mesh->getProperties().addProperty(prop_name_3, properties_3,
		MeshLib::MeshItemType::Cell);

	// fetch the vector in order to compare the content
	boost::optional<
		MeshLib::PropertyVector<Eigen::Matrix<double,3,3,Eigen::RowMajor>>*
	> properties_3_cpy(
		mesh->getProperties().getProperty<Eigen::Matrix<double,3,3,Eigen::RowMajor>>(
			prop_name_3, MeshLib::MeshItemType::Cell
		)
	);
	ASSERT_FALSE(!properties_3_cpy);

	// compare the values/matrices
	auto it_cpy=(*properties_3_cpy)->begin();
	for (auto it=properties_3->begin(); it != properties_3->end();
		it++, it_cpy++) {
		for (int r(0); r<it->rows(); r++) {
			for (int c(0); c<it->cols(); c++) {
				ASSERT_EQ((*it)(r,c), (*it_cpy)(r,c));
			}
		}
	}
#endif
}

