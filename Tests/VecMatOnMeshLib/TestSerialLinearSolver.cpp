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

#include <vector>
#include <memory>
#include <cmath>

#include <gtest/gtest.h>

#include "MathLib/MathTools.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/LinAlg/Dense/DenseTools.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"

#include "VecMatOnMeshLib/VecMeshItems/VectorComposition.h"
#include "VecMatOnMeshLib/VecMeshItems/MeshItem.h"
#include "VecMatOnMeshLib/MeshItemWiseTask/LinearSystemAssembler.h"
#include "VecMatOnMeshLib/Serial/SerialVecMatOnMesh.h"
#include "VecMatOnMeshLib/Serial/ForEachMeshItem.h"

#include "../TestTools.h"
#include "SteadyDiffusion2DExample1.h"


TEST(VecMatOnMeshLib, SerialLinearSolver)
{
    // example
    SteadyDiffusion2DExample1 ex1;
    typedef SteadyDiffusion2DExample1::LocalAssembler MyLocalAssembler;
    MyLocalAssembler local_2d_diff;

    //--------------------------------------------------------------------------
    // Choose implementation type
    //--------------------------------------------------------------------------
    typedef VecMatOnMeshLib::SerialVecMatOnMesh TVecMatOnMesh;
    typedef TVecMatOnMesh::VectorType TVec;
    typedef TVecMatOnMesh::MatrixType TMat;
    TVecMatOnMesh vecMatOnMesh;

    //--------------------------------------------------------------------------
    // Prepare mesh items where data are assigned
    //--------------------------------------------------------------------------
    const VecMatOnMeshLib::MeshItems mesh_items_all_nodes(*ex1.msh, ex1.msh->getNodes());

    //--------------------------------------------------------------------------
    // Allocate a coefficient matrix, RHS and solution vectors
    //--------------------------------------------------------------------------
    // define a mesh item composition in a vector
    std::vector<VecMatOnMeshLib::ComponentDistribution*> vec_comp_dis;
    vec_comp_dis.push_back(new VecMatOnMeshLib::ComponentDistribution(&mesh_items_all_nodes));
    VecMatOnMeshLib::VectorComposition vec1_composition(vec_comp_dis, VecMatOnMeshLib::OrderingType::BY_COMPONENT_TYPE);

    // allocate a vector and matrix
    std::unique_ptr<TMat> A(vecMatOnMesh.createMatrix(vec1_composition));
    std::unique_ptr<TVec> rhs(vecMatOnMesh.createVector(vec1_composition));
    std::unique_ptr<TVec> x(vecMatOnMesh.createVector(vec1_composition));

    //--------------------------------------------------------------------------
    // Construct a linear system
    //--------------------------------------------------------------------------
    // create a mapping table from element nodes to entries in the linear system
    auto &all_eles = ex1.msh->getElements();
    std::vector<std::vector<std::size_t> > map_ele_nodes2vec_entries(all_eles.size());
    for (std::size_t i=0; i<map_ele_nodes2vec_entries.size(); i++) {
        auto* e = all_eles[i];
        std::vector<VecMatOnMeshLib::MeshItem> vec_items;
        for (std::size_t j=0; j<e->getNNodes(); j++)
            vec_items.push_back(VecMatOnMeshLib::MeshItem(ex1.msh->getID(), VecMatOnMeshLib::MeshItemType::Node, e->getNode(j)->getID()));
        map_ele_nodes2vec_entries[i] = vec1_composition.getDataIDList(vec_items, VecMatOnMeshLib::OrderingType::BY_COMPONENT_TYPE);
        //std::cout << i << ": "; std::for_each(map_ele_nodes2vec_entries[i].begin(), map_ele_nodes2vec_entries[i].end(), [](std::size_t id){std::cout << id << " ";}); std::cout << "\n";
    }

    // create a local assembler
    typedef VecMatOnMeshLib::LinearSystemAssembler<TMat, TVec, MeshLib::Element, MyLocalAssembler> LocalAssembler;
    LocalAssembler assembler(*A.get(), *rhs.get(), local_2d_diff, map_ele_nodes2vec_entries);

    // do assembly
    TVecMatOnMesh::ForEachType<MeshLib::Element, LocalAssembler> vec1_global_assembly;
    vec1_global_assembly(ex1.msh->getElements(), assembler);
    //std::cout << "A=\n"; A->write(std::cout); std::cout << "rhs=\n"; rhs->write(std::cout);

    // apply Dirichlet BC
    MathLib::applyKnownSolution(*A, *rhs, ex1.vec_DirichletBC_id, ex1.vec_DirichletBC_value);
    //std::cout << "A=\n"; A->write(std::cout); std::cout << "rhs=\n"; rhs->write(std::cout);

    //--------------------------------------------------------------------------
    // solve x=A^-1 rhs
    //--------------------------------------------------------------------------
    MathLib::GaussAlgorithm<TMat, TVec> ls(*A);
    ls.solve(*rhs, *x);

    double *px = &(*x)[0];
    ASSERT_DOUBLE_ARRAY_EQ(&ex1.exact_solutions[0], px, 9, 1.e-5);
}
