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

#include "NumLib/Assembler/VectorMatrixAssembler.h"

#include "MathLib/LinAlg/ApplyKnownSolution.h"
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "MathLib/LinAlg/MatrixSpecifications.h"
#include "MathLib/LinAlg/MatrixVectorTraits.h"
#include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
#include "MathLib/LinAlg/FinalizeMatrixAssembly.h"
#include "MathLib/MathTools.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Location.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

#include "NumLib/NumericsConfig.h"

#include "../TestTools.h"
#include "SteadyDiffusion2DExample1.h"

TEST(NumLibSerialLinearSolver, Steady2DdiffusionQuadElem)
{
    // example
    using Example = SteadyDiffusion2DExample1<GlobalIndexType>;
    Example ex1;

    //--------------------------------------------------------------------------
    // Prepare mesh items where data are assigned
    //--------------------------------------------------------------------------
    const MeshLib::MeshSubset mesh_items_all_nodes(*ex1.msh,
                                                          &ex1.msh->getNodes());

    //--------------------------------------------------------------------------
    // Allocate a coefficient matrix, RHS and solution vectors
    //--------------------------------------------------------------------------
    // define a mesh item composition in a vector
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> vec_comp_dis;
    vec_comp_dis.emplace_back(new MeshLib::MeshSubsets{&mesh_items_all_nodes});
    NumLib::LocalToGlobalIndexMap local_to_global_index_map(
        std::move(vec_comp_dis), NumLib::ComponentOrder::BY_COMPONENT);

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
    auto x = MathLib::MatrixVectorTraits<GlobalVector>::newInstance(ms);
    // TODO no setZero() for rhs, x?

    using LocalAssembler = Example::LocalAssemblerData;
    // Initializer of the local assembler data.
    std::vector<LocalAssembler*> local_assembler_data;
    local_assembler_data.resize(ex1.msh->getNumberOfElements());

    auto local_asm_builder =
        [&](std::size_t const id,
            MeshLib::Element const& item,
            LocalAssembler*& item_data)
    {
        assert(local_to_global_index_map.size() > id);

        auto const num_local_dof = local_to_global_index_map.getNumElementDOF(id);

        Example::initializeLocalData(
                    item, item_data, num_local_dof, ex1);
    };

    // Call global initializer for each mesh element.
    GlobalExecutor::transformDereferenced(
            local_asm_builder,
            ex1.msh->getElements(),
            local_assembler_data);

    // TODO in the future use simpler NumLib::ODESystemTag
    // Local and global assemblers.
    typedef NumLib::VectorMatrixAssembler<
            LocalAssembler,
            NumLib::ODESystemTag::FirstOrderImplicitQuasilinear> GlobalAssembler;

    GlobalAssembler assembler(local_to_global_index_map);

    // Call global assembler for each mesh element.
    auto M_dummy = MathLib::MatrixVectorTraits<GlobalMatrix>::newInstance(ms);
    A->setZero();
    auto const t = 0.0;
    GlobalExecutor::executeMemberDereferenced(
                assembler, &GlobalAssembler::assemble,
                local_assembler_data, t, *x, *M_dummy, *A, *rhs);

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
    for (std::size_t i = 0; i < x->size(); ++i)
        solution[i] = (*x)[i];

    ASSERT_ARRAY_NEAR(&ex1.exact_solutions[0], &solution[0], ex1.dim_eqs, 1.e-5);

    for (auto p : local_assembler_data)
        delete p;
}
