/**
 * \author Norihiro Watanabe
 * \date   2013-04-18
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef STEADYDIFFUSION2DEXAMPLE1_H_
#define STEADYDIFFUSION2DEXAMPLE1_H_

#include <cmath>
#include <vector>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Node.h"

struct SteadyDiffusion2DExample1
{
	class LocalAssembler
	{
	public:
		LocalAssembler()
			: _m(4,4)
		{
			_m(0,0) = 4.0; _m(0,1) = -1.0; _m(0,2) = -2.0; _m(0,3) = -1.0;
			_m(1,1) = 4.0; _m(1,2) = -1.0; _m(1,3) = -2.0;
			_m(2,2) = 4.0; _m(2,3) = -1.0;
			_m(3,3) = 4.0;

			// copy upper triangle to lower to make symmetric
			for (std::size_t i = 0; i < 4; i++)
				for (std::size_t j = 0; j < i; j++)
					_m(i,j) = _m(j,i);

			//_m *= 1.e-11/6.0;
			for (std::size_t i = 0; i < 4; i++)
				for (std::size_t j = 0; j < 4; j++)
					_m(i,j) *= 1.e-11 / 6.0;
		}

		void operator()(const MeshLib::Element & /*e*/,
		                MathLib::DenseMatrix<double> &localA,
		                MathLib::DenseVector<double> & /*rhs*/)
		{
			for (std::size_t i = 0; i < 4; i++)
				for (std::size_t j = 0; j < 4; j++)
					localA(i,j) = _m(i,j);
		}

	private:
		MathLib::DenseMatrix<double> _m;
	};

	SteadyDiffusion2DExample1()
	{
		msh = MeshLib::MeshGenerator::generateRegularQuadMesh(2.0, mesh_subdivs);
		for (auto itr=msh->getNodes().cbegin(); itr!=msh->getNodes().cend(); ++itr) {
			auto* node = *itr;
			vec_nodeIDs.push_back(node->getID());
		}
		vec_DirichletBC_id.resize(2 * mesh_stride);
		for (std::size_t i = 0; i < mesh_stride; i++)
		{
			// left side
			vec_DirichletBC_id[i] = i * mesh_stride;

			// right side
			vec_DirichletBC_id[mesh_stride + i] = i * mesh_stride + mesh_subdivs;
		}

		vec_DirichletBC_value.resize(2 * mesh_stride);
		std::fill_n(vec_DirichletBC_value.begin(), mesh_stride, 0);
		std::fill_n(vec_DirichletBC_value.begin() + mesh_stride, mesh_stride, 1);
		exact_solutions.resize(dim_eqs);
		for (std::size_t i = 0; i < mesh_stride; i++)
			for (std::size_t j = 0; j < mesh_stride; j++)
				exact_solutions[i*mesh_stride + j] = j * 1./mesh_subdivs;
	}

	~SteadyDiffusion2DExample1()
	{
		delete msh;
	}

	static const std::size_t mesh_subdivs = 10;
	static const std::size_t mesh_stride = mesh_subdivs + 1;
	static const std::size_t dim_eqs = mesh_stride * mesh_stride;
	MeshLib::Mesh* msh;
	std::vector<std::size_t> vec_DirichletBC_id;
	std::vector<double> vec_DirichletBC_value;
	std::vector<double> exact_solutions;
	std::vector<std::size_t> vec_nodeIDs;
};

#endif //STEADYDIFFUSION2DEXAMPLE1_H_
