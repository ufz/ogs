/**
 * \author Norihiro Watanabe
 * \date   2013-04-18
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "AssemblerLib/VectorMatrixAssembler.h"
#include "AssemblerLib/LocalAssemblerBuilder.h"
#include "AssemblerLib/SerialDenseSetup.h"


#include "MathLib/LinAlg/Dense/DenseTools.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/MathTools.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Location.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

#include "../TestTools.h"
#include "SteadyDiffusion2DExample1.h"

TEST(AssemblerLibSerialLinearSolver, Steady2DdiffusionQuadElem)
{
    // example
    SteadyDiffusion2DExample1 ex1;

    //--------------------------------------------------------------------------
    // Choose implementation type
    //--------------------------------------------------------------------------
    typedef AssemblerLib::SerialDenseSetup GlobalSetup;
    const GlobalSetup globalSetup;

    //--------------------------------------------------------------------------
    // Prepare mesh items where data are assigned
    //--------------------------------------------------------------------------
    const MeshLib::MeshSubset mesh_items_all_nodes(*ex1.msh,
                                                           ex1.msh->getNodes());

    //--------------------------------------------------------------------------
    // Allocate a coefficient matrix, RHS and solution vectors
    //--------------------------------------------------------------------------
    // define a mesh item composition in a vector
    std::vector<MeshLib::MeshSubsets*> vec_comp_dis;
    vec_comp_dis.push_back(
        new MeshLib::MeshSubsets(&mesh_items_all_nodes));
    AssemblerLib::LocalToGlobalIndexMap local_to_global_index_map(
            vec_comp_dis);

    //--------------------------------------------------------------------------
    // Construct a linear system
    //--------------------------------------------------------------------------
    // allocate a vector and matrix
    typedef GlobalSetup::VectorType GlobalVector;
    typedef GlobalSetup::MatrixType GlobalMatrix;
    std::unique_ptr<GlobalMatrix> A(globalSetup.createMatrix(local_to_global_index_map.dofSize()));
    A->setZero();
    std::unique_ptr<GlobalVector> rhs(globalSetup.createVector(local_to_global_index_map.dofSize()));
    std::unique_ptr<GlobalVector> x(globalSetup.createVector(local_to_global_index_map.dofSize()));

    // Initializer of the local assembler data.
    std::vector<SteadyDiffusion2DExample1::LocalAssemblerData*>
        local_assembler_data;
    local_assembler_data.resize(ex1.msh->getNElements());

    typedef AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            void (const MeshLib::Element &,
                    SteadyDiffusion2DExample1::LocalAssemblerData *&,
                    SteadyDiffusion2DExample1 const&)
                > LocalAssemblerBuilder;

    LocalAssemblerBuilder local_asm_builder(
        SteadyDiffusion2DExample1::initializeLocalData);

    // Call global initializer for each mesh element.
    globalSetup.execute(
            local_asm_builder,
            ex1.msh->getElements(),
            local_assembler_data,
            ex1);

    // Local and global assemblers.
    typedef AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector> GlobalAssembler;

    GlobalAssembler assembler(*A.get(), *rhs.get(), local_to_global_index_map);

    // Call global assembler for each mesh element.
    globalSetup.execute(assembler, local_assembler_data);

    //std::cout << "A=\n";
    //A->write(std::cout);
    //std::cout << "rhs=\n";
    //rhs->write(std::cout);

    // apply Dirichlet BC
    MathLib::applyKnownSolution(*A, *rhs, ex1.vec_DirichletBC_id,
                                ex1.vec_DirichletBC_value);
    //std::cout << "A=\n";
    //A->write(std::cout);
    //std::cout << "rhs=\n";
    //rhs->write(std::cout);

    MathLib::finalizeMatrixAssembly(*A);
    //--------------------------------------------------------------------------
    // solve x=A^-1 rhs
    //--------------------------------------------------------------------------
    MathLib::GaussAlgorithm<GlobalMatrix, GlobalVector> ls(*A);
    ls.solve(*rhs, *x);

    double* px = &(*x)[0];
    ASSERT_ARRAY_NEAR(&ex1.exact_solutions[0], px, ex1.dim_eqs, 1.e-5);

    std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
        [](MeshLib::MeshSubsets * p) { delete p; return true; });
}
