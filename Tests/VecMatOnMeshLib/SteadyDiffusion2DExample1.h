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

#ifndef STEADYDIFFUSION2DEXAMPLE1_H_
#define STEADYDIFFUSION2DEXAMPLE1_H_

#include <cmath>
#include <vector>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/DenseVector.h"

#include "MeshLib/Elements/Edge.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
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
		msh = MeshLib::MeshGenerator::generateRegularQuadMesh(2.0, 2);
		for (auto* node : msh->getNodes())
			vec_nodeIDs.push_back(node->getID());
		std::size_t int_dirichlet_bc_id[] = {2,5,8,0,3,6};
		vec_DirichletBC_id.assign(int_dirichlet_bc_id,
		                          int_dirichlet_bc_id + 6);
		vec_DirichletBC_value.resize(6);
		std::fill_n(vec_DirichletBC_value.begin(), 3, 0);
		std::fill_n(vec_DirichletBC_value.begin() + 3, 3, 1);
		exact_solutions.resize(9);
		for (std::size_t i = 0; i < 9; i++)
		{
			if (i % 3 == 0) exact_solutions[i] = 1.0;
			if (i % 3 == 1) exact_solutions[i] = 0.5;
			if (i % 3 == 2) exact_solutions[i] = 0.;
		}
	}

	~SteadyDiffusion2DExample1()
	{
		delete msh;
	}

	static const std::size_t dim_eqs = 9;
	MeshLib::Mesh* msh;
	std::vector<std::size_t> vec_DirichletBC_id;
	std::vector<double> vec_DirichletBC_value;
	std::vector<double> exact_solutions;
	std::vector<std::size_t> vec_nodeIDs;
}; //SteadyDiffusion2DExample1

#endif //STEADYDIFFUSION2DEXAMPLE1_H_
