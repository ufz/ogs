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
#include "AssemblerLib/MeshComponentMap.h"
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
    AssemblerLib::MeshComponentMap vec1_composition(
        vec_comp_dis, AssemblerLib::ComponentOrder::BY_COMPONENT);

    // allocate a vector and matrix
    typedef GlobalSetup::VectorType GlobalVector;
    typedef GlobalSetup::MatrixType GlobalMatrix;
    std::unique_ptr<GlobalMatrix> A(globalSetup.createMatrix(vec1_composition));
    A->setZero();
    std::unique_ptr<GlobalVector> rhs(globalSetup.createVector(vec1_composition));
    std::unique_ptr<GlobalVector> x(globalSetup.createVector(vec1_composition));

    //--------------------------------------------------------------------------
    // Construct a linear system
    //--------------------------------------------------------------------------
    // create a mapping table from element nodes to entries in the linear system
    auto &all_eles = ex1.msh->getElements();
    std::vector<std::vector<std::size_t> > map_ele_nodes2vec_entries;
    map_ele_nodes2vec_entries.reserve(all_eles.size());
    for (auto e = all_eles.cbegin(); e != all_eles.cend(); ++e)
    {
        std::size_t const nnodes = (*e)->getNNodes();
        std::size_t const mesh_id = ex1.msh->getID();
        std::vector<MeshLib::Location> vec_items;
        vec_items.reserve(nnodes);
        for (std::size_t j = 0; j < nnodes; j++)
            vec_items.emplace_back(
                mesh_id,
                MeshLib::MeshItemType::Node,
                (*e)->getNode(j)->getID());

        map_ele_nodes2vec_entries.push_back(
            vec1_composition.getGlobalIndices
                <AssemblerLib::ComponentOrder::BY_COMPONENT>(vec_items));
    }

    // Initializer of the local assembler data.
    using LocalMatrixType = MathLib::DenseMatrix<double>;
    using LocalVectorType = MathLib::DenseVector<double>;

    std::vector<LocalAssemblerData<LocalMatrixType, LocalVectorType>*>
        local_assembler_data;
    local_assembler_data.resize(ex1.msh->getNElements());

    LocalDataInitializer<LocalMatrixType, LocalVectorType> initializer;

    typedef AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            LocalDataInitializer<LocalMatrixType, LocalVectorType>
                > GlobalInitializer;

    GlobalInitializer global_initializer(initializer);

    // Call global initializer for each mesh element.
    globalSetup.execute(
            global_initializer,
            ex1.msh->getElements(),
            local_assembler_data,
            ex1._localA,
            ex1._localRhs);

    // Local and global assemblers.
    typedef AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector> GlobalAssembler;

    GlobalAssembler assembler(*A.get(), *rhs.get(),
        AssemblerLib::LocalToGlobalIndexMap(map_ele_nodes2vec_entries));

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
