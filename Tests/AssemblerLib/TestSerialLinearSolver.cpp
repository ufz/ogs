/**
 * \author Norihiro Watanabe
 * \date   2013-04-18
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/MathTools.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Location.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

#include "ProcessLib/NumericsConfig.h"

#include "../TestTools.h"
#include "SteadyDiffusion2DExample1.h"

TEST(AssemblerLibSerialLinearSolver, Steady2DdiffusionQuadElem)
{
    // example
    using Example = SteadyDiffusion2DExample1<GlobalIndexType>;
    Example ex1;

    //--------------------------------------------------------------------------
    // Choose implementation type
    //--------------------------------------------------------------------------
    using GlobalSetup = GlobalSetupType;    // defined in numerics config
    const GlobalSetup globalSetup;

    //--------------------------------------------------------------------------
    // Prepare mesh items where data are assigned
    //--------------------------------------------------------------------------
    const MeshLib::MeshSubset mesh_items_all_nodes(*ex1.msh,
                                                          &ex1.msh->getNodes());

    //--------------------------------------------------------------------------
    // Allocate a coefficient matrix, RHS and solution vectors
    //--------------------------------------------------------------------------
    // define a mesh item composition in a vector
    std::vector<MeshLib::MeshSubsets*> vec_comp_dis;
    vec_comp_dis.push_back(
        new MeshLib::MeshSubsets(&mesh_items_all_nodes));
    AssemblerLib::LocalToGlobalIndexMap local_to_global_index_map(
            vec_comp_dis, AssemblerLib::ComponentOrder::BY_COMPONENT);

    //--------------------------------------------------------------------------
    // Construct a linear system
    //--------------------------------------------------------------------------
    // allocate a vector and matrix
    typedef GlobalSetup::VectorType GlobalVector;
    typedef GlobalSetup::MatrixType GlobalMatrix;
    auto A   = globalSetup.createMatrix(local_to_global_index_map.dofSize());
    A->setZero();
    auto rhs = globalSetup.createVector(local_to_global_index_map.dofSize());
    auto x   = globalSetup.createVector(local_to_global_index_map.dofSize());
    // TODO no setZero() for rhs, x?

    // Initializer of the local assembler data.
    std::vector<Example::LocalAssemblerData<
        GlobalMatrix, GlobalVector>*> local_assembler_data;
    local_assembler_data.resize(ex1.msh->getNElements());

    typedef AssemblerLib::LocalAssemblerBuilder<
            MeshLib::Element,
            void (const MeshLib::Element &,
                    Example::LocalAssemblerData<
                        GlobalMatrix, GlobalVector> *&,
                    std::size_t const local_matrix_size,
                    Example const&)
                > LocalAssemblerBuilder;

    LocalAssemblerBuilder local_asm_builder(
        Example::initializeLocalData<GlobalMatrix, GlobalVector>,
        local_to_global_index_map);

    // Call global initializer for each mesh element.
    globalSetup.transform(
            local_asm_builder,
            ex1.msh->getElements(),
            local_assembler_data,
            ex1);

    // TODO in the future use simpler NumLib::ODESystemTag
    // Local and global assemblers.
    typedef AssemblerLib::VectorMatrixAssembler<
            GlobalMatrix, GlobalVector,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear> GlobalAssembler;

    GlobalAssembler assembler(local_to_global_index_map);

    // Call global assembler for each mesh element.
    auto M_dummy = globalSetup.createMatrix(local_to_global_index_map.dofSize());
    A->setZero();
    auto const t = 0.0;
    globalSetup.execute(assembler, local_assembler_data, t, *x, *M_dummy, *A, *rhs);

    //std::cout << "A=\n";
    //A->write(std::cout);
    //std::cout << "rhs=\n";
    //rhs->write(std::cout);

    // apply Dirichlet BC
    MathLib::applyKnownSolution(*A, *rhs, *x, ex1.vec_DirichletBC_id,
                                ex1.vec_DirichletBC_value);
    //std::cout << "A=\n";
    //A->write(std::cout);
    //std::cout << "rhs=\n";
    //rhs->write(std::cout);

    MathLib::finalizeMatrixAssembly(*A);
    //--------------------------------------------------------------------------
    // solve x=A^-1 rhs
    //--------------------------------------------------------------------------
    boost::property_tree::ptree t_root;
    {
        boost::property_tree::ptree t_solver;
        t_solver.put("solver_type", "CG");
        t_solver.put("precon_type", "NONE");
        t_solver.put("error_tolerance", 1e-16);
        t_solver.put("max_iteration_step", 1000);
        t_root.put_child("eigen", t_solver);
    }
    t_root.put("lis", "-i cg -p none -tol 1e-16 -maxiter 1000");
    BaseLib::ConfigTree conf(t_root, "",
                             BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    GlobalSetup::LinearSolver ls("solver_name", &conf);
    ls.solve(*A, *rhs, *x);

    // copy solution to double vector
    std::vector<double> solution(x->size());
    for (std::size_t i = 0; i < x->size(); ++i)
        solution[i] = (*x)[i];

    ASSERT_ARRAY_NEAR(&ex1.exact_solutions[0], &solution[0], ex1.dim_eqs, 1.e-5);

    std::remove_if(vec_comp_dis.begin(), vec_comp_dis.end(),
        [](MeshLib::MeshSubsets * p) { delete p; return true; });

    for (auto p : local_assembler_data)
        delete p;
}
