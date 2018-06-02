/**
 * \author Norihiro Watanabe
 * \date   2013-04-18
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/MatrixSpecifications.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/MathTools.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Location.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/Node.h"

#include "NumLib/NumericsConfig.h"
#include "NumLib/DOF/DOFTableUtil.h"

#include "Tests/TestTools.h"
#include "SteadyDiffusion2DExample1.h"

TEST(NumLibSerialLinearSolver, Steady2DdiffusionQuadElem)
{
    // example
    using Example = SteadyDiffusion2DExample1<GlobalIndexType>;
    Example ex1;

    //--------------------------------------------------------------------------
    // Prepare mesh items where data are assigned
    //--------------------------------------------------------------------------
    MeshLib::MeshSubset const mesh_items_all_nodes{*ex1.msh,
                                                   ex1.msh->getNodes()};

    //--------------------------------------------------------------------------
    // Allocate a coefficient matrix, RHS and solution vectors
    //--------------------------------------------------------------------------
    // define a mesh item composition in a vector
    NumLib::LocalToGlobalIndexMap local_to_global_index_map(
        {mesh_items_all_nodes}, NumLib::ComponentOrder::BY_COMPONENT);

    //--------------------------------------------------------------------------
    // Construct a linear system
    //--------------------------------------------------------------------------
    // allocate a vector and matrix
    MathLib::MatrixSpecifications ms{local_to_global_index_map.dofSizeWithoutGhosts(),
        local_to_global_index_map.dofSizeWithoutGhosts(),
        &local_to_global_index_map.getGhostIndices(),
        nullptr};
    auto A = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(ms);
    A->setZero();
    auto rhs = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(ms);
    rhs->setZero();
    auto x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(ms);
    x->setZero();

    using LocalAssembler = Example::LocalAssemblerData;
    // Initializer of the local assembler data.
    std::vector<LocalAssembler*> local_assemblers;
    local_assemblers.resize(ex1.msh->getNumberOfElements());

    auto local_asm_builder =
        [&](std::size_t const id,
            MeshLib::Element const& item,
            LocalAssembler*& item_data)
    {
        assert(local_to_global_index_map.size() > id);

        auto const num_local_dof = local_to_global_index_map.getNumberOfElementDOF(id);

        Example::initializeLocalData(
                    item, item_data, num_local_dof, ex1);
    };

    // Call global initializer for each mesh element.
    GlobalExecutor::transformDereferenced(
            local_asm_builder,
            ex1.msh->getElements(),
            local_assemblers);

    // Call global assembler for each mesh element.
    auto M_dummy = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(ms);
    A->setZero();
    auto const t = 0.0;
    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssembler::assemble, local_assemblers, local_to_global_index_map,
        t, *x, *M_dummy, *A, *rhs);

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

    GlobalLinearSolver ls("solver_name", &conf);
    ls.solve(*A, *rhs, *x);

    // copy solution to double vector
    std::vector<double> solution(x->size());
    for (GlobalIndexType i = 0; i < x->size(); ++i)
        solution[i] = (*x)[i];

    ASSERT_ARRAY_NEAR(&ex1.exact_solutions[0], &solution[0], ex1.dim_eqs, 1.e-5);

    for (auto p : local_assemblers)
        delete p;
}
