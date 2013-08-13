/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-18
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseTools.h"
#include "MathLib/MathTools.h"

#include "MeshLib/Elements/Edge.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Location.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"

#include "VecMatOnMeshLib/MeshItemWiseTask/LocalToGlobalIndexMap.h"
#include "VecMatOnMeshLib/MeshItemWiseTask/MatrixAssembler.h"
#include "VecMatOnMeshLib/MeshItemWiseTask/VectorAssembler.h"
#include "VecMatOnMeshLib/Serial/SerialVectorMatrixBuilder.h"
#include "VecMatOnMeshLib/VecMeshItems/MeshComponentMap.h"

#include "../TestTools.h"

namespace
{
/// Local assembler for creating a vector of x coords of nodes
class LocalVecAssemblerExtractNodeX
{
public:
	void operator()(const MeshLib::Node &node,
	                MathLib::DenseVector<double> &local_vec)
	{
		local_vec[0] = node[0]; // set x
	}
};

/// Local assembler for creating a matrix having distances between element nodes
class LocalMatAssemblerNodeDistance
{
public:
	template <class T_ELE>
	void operator()(const T_ELE &element,
	                MathLib::DenseMatrix<double> &local_mat)
	{
		assert(local_mat.getNRows() >= element.getNNodes()
		    && local_mat.getNCols() >= element.getNNodes());

		const std::size_t nnodes = element.getNNodes();
		for (std::size_t i = 0; i < nnodes; i++)
			for (std::size_t j = 0; j < nnodes; j++)
				local_mat(i, j) = std::sqrt(
				    MathLib::sqrDist(element.getNode(i), element.getNode(j)));

	}
};

} // end anonymous namespace

TEST(VecMatOnMeshLib, SerialVectorMatrixBuilder)
{
	// This test case checks:
	// - construct a vector "v" having x coords. of nodes in the left-half domain
	// - construct a linear operator matrix "m" for "v" (m has node distances)
	// - calculate r = m*v

	//--------------------------------------------------------------------------
	// Choose implementation type
	//--------------------------------------------------------------------------
	typedef VecMatOnMeshLib::SerialDenseVectorMatrixBuilder SerialBuilder;
	typedef SerialBuilder::VectorType TVec;
	typedef SerialBuilder::MatrixType TMat;
	SerialBuilder vecMatOnMesh;

	//--------------------------------------------------------------------------
	// Prepare a mesh having line elements
	//--------------------------------------------------------------------------
	// create a line mesh with 10 elements and 11 nodes (dx=0.1)
	std::unique_ptr<MeshLib::Mesh> msh(
	    MeshLib::MeshGenerator::generateLineMesh(1.0, 10));

	//--------------------------------------------------------------------------
	// Prepare mesh items (nodes and elements) where data are assigned
	//--------------------------------------------------------------------------
	// assign data to nodes in left-half of the domain
	// (x<=0.5, corresponding to first six nodes)
	std::vector<MeshLib::Node*> vec_selected_nodes(6);
	{
		auto &vec_all_nodes = msh->getNodes();
		for (std::size_t i = 0; i < vec_selected_nodes.size(); i++)
			vec_selected_nodes[i] = vec_all_nodes[i];
	}
	const MeshLib::MeshSubset mesh_items_left_nodes(*msh,
	                                                        vec_selected_nodes);

	//extract elements having those nodes (corresponding to first five elements)
	std::vector<MeshLib::Element*> vec_selected_eles(5);
	auto &vec_all_eles = msh->getElements();
	for (std::size_t i = 0; i < vec_selected_eles.size(); i++)
		vec_selected_eles[i] = vec_all_eles[i];

	//--------------------------------------------------------------------------
	// Allocate a vector and a linear operator (i.e. matrix)
	//--------------------------------------------------------------------------
	// define a mesh item composition in a vector
	std::vector<MeshLib::MeshSubsets*> vec_comp_dis;
	vec_comp_dis.push_back(
	    new MeshLib::MeshSubsets(&mesh_items_left_nodes));

	VecMatOnMeshLib::MeshComponentMap vec1_composition(
	    vec_comp_dis, VecMatOnMeshLib::ComponentOrder::BY_COMPONENT);
	//vec1_composition.print();

	// allocate a vector and matrix
	std::unique_ptr<TVec> vec_left_nodes_x_coord(
	    vecMatOnMesh.createVector(vec1_composition));
	std::unique_ptr<TMat> mat(vecMatOnMesh.createMatrix(vec1_composition));
	mat->setZero();

	//--------------------------------------------------------------------------
	// Construct the vector from the selected nodes
	//--------------------------------------------------------------------------
	// prepare a mapping table from DoFs to positions in the vector
	// node id, comp id <-> vector entry id
	std::vector<std::vector<std::size_t> > map_node2vec_entry(
	        mesh_items_left_nodes.getNNodes());
	for (std::size_t i = 0; i < map_node2vec_entry.size(); i++)
		map_node2vec_entry[i] = vec1_composition.getDataIDList(
	        MeshLib::Location(msh->getID(),
	                                  MeshLib::MeshItemType::Node,
		                              mesh_items_left_nodes.getNodeID(i)));

	// create a vector assembler
	LocalVecAssemblerExtractNodeX extractX;
	typedef VecMatOnMeshLib::VectorAssembler<
	    TVec, MeshLib::Node, LocalVecAssemblerExtractNodeX> SetNodeXToVec;
	SetNodeXToVec vec1Assembler(
	    *vec_left_nodes_x_coord.get(), extractX, map_node2vec_entry);

	// do assembly for each selected node
	SerialBuilder::forEachMeshItem(vec_selected_nodes, vec1Assembler);

	ASSERT_EQ(0.0, (*vec_left_nodes_x_coord)[0]);
	ASSERT_EQ(0.5, (*vec_left_nodes_x_coord)[5]);

	//--------------------------------------------------------------------------
	// Construct the matrix from the selected elements
	//--------------------------------------------------------------------------
	// prepare a mapping table from DoFs to positions in the matrix
	// element id, node id, comp id <-> vector entry id
	std::vector<std::vector<std::size_t> > mat_row_column_positions;
	{
		mat_row_column_positions.reserve(vec_selected_eles.size());
		for (std::size_t i = 0; i < vec_selected_eles.size(); i++)
		{
			auto* e = vec_selected_eles[i];
			std::vector<MeshLib::Location> vec_items;
			for (std::size_t j = 0; j < e->getNNodes(); j++)
				vec_items.push_back(
					MeshLib::Location(msh->getID(),
					MeshLib::MeshItemType::Node,
					e->getNode(j)->getID()));

			mat_row_column_positions.push_back(vec1_composition.getDataIDList
				<VecMatOnMeshLib::ComponentOrder::BY_COMPONENT>(vec_items));
			//std::cout << i << ": ";
			//std::copy(vec_data_pos[i].begin(), vec_data_pos[i].end(),
			//    std::ostream_iterator<std::size_t>(std::cout, " "));
			//    std::cout << "\n";
		}
	}
	VecMatOnMeshLib::LocalToGlobalIndexMap mat_data_pos(mat_row_column_positions);


	// create a matrix assembler
	LocalMatAssemblerNodeDistance nodeDist;
	typedef VecMatOnMeshLib::MatrixAssembler<
	    TMat, MeshLib::Element, LocalMatAssemblerNodeDistance> SetNodeDistXToMat;
	SetNodeDistXToMat matAssembler(*mat.get(), nodeDist, mat_data_pos);

	// do assembly for each selected element
	SerialBuilder::forEachMeshItem(vec_selected_eles, matAssembler);

	ASSERT_NEAR(0.0, (*mat)(0,0), 1e-6);
	ASSERT_NEAR(0.1, (*mat)(0,1), 1e-6);
	ASSERT_NEAR(0.1, (*mat)(1,0), 1e-6);
	ASSERT_NEAR(0.0, (*mat)(1,1), 1e-6);
	ASSERT_NEAR(0.1, (*mat)(1,2), 1e-6);
	ASSERT_NEAR(0.1, (*mat)(4,3), 1e-6);
	ASSERT_NEAR(0., (*mat)(4,4), 1e-6);

	//--------------------------------------------------------------------------
	// an example of using the constructed matrix and vector
	// y = A * x
	//--------------------------------------------------------------------------
	std::unique_ptr<TVec> vec_out(vecMatOnMesh.createVector(vec1_composition));
	mat->matvec(*vec_left_nodes_x_coord, *vec_out);

	ASSERT_NEAR(0.01, (*vec_out)[0], 1e-6);
	ASSERT_NEAR(0.02, (*vec_out)[1], 1e-6);
	ASSERT_NEAR(0.04, (*vec_out)[2], 1e-6);

	std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
		[](MeshLib::MeshSubsets * p) { delete p; return true; });
}
