#pragma once

#if 0

// Moved here from TESProcess.cpp
template<typename GlobalSetup>
void TESProcess<GlobalSetup>::
assembleConcreteProcess(
        const double t, GlobalVector const& x,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble TESProcess.");

    unsigned num_try = 0;

    INFO("-> TES process try number %u in current picard iteration", num_try);
    _assembly_params.number_of_try_of_iteration = num_try;

    // Call global assembler for each local assembly item.
    _global_setup.execute(*BP::_global_assembler, _local_assemblers,
                          t, x, M, K, b);

#if 0 && defined(OGS_USE_EIGEN) && ! defined(OGS_USE_EIGENLIS)
    // TODO put that somewhere
    MathLib::scaleDiagonal(*_A, *_rhs);
#endif

#if 0 && defined(OGS_USE_EIGENLIS)
    // TODO put that somewhere

    // scaling
    typename GlobalMatrix::RawMatrixType AT = _A->getRawMatrix().transpose();

    for (unsigned dof = 0; dof < NODAL_DOF; ++dof)
    {
        auto const& trafo = (dof == 0) ? _assembly_params.trafo_p
                          : (dof == 1) ? _assembly_params.trafo_T
                                       : _assembly_params.trafo_x;

        for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
        {
            MeshLib::Location loc(BP::_mesh.getID(), MeshLib::MeshItemType::Node, i);
            auto const idx = BP::_local_to_global_index_map->getGlobalIndex(loc, dof);

            AT.row(idx) *= trafo.dxdy(0);
            x_curr[idx] /= trafo.dxdy(0);
        }
    }

    _A->getRawMatrix() = AT.transpose();
#endif

#ifndef NDEBUG
    if (_total_iteration == 0 && num_try == 0 && _output_global_matrix)
    {
        MathLib::BLAS::finalizeAssembly(M);
        MathLib::BLAS::finalizeAssembly(K);
        MathLib::BLAS::finalizeAssembly(b);

        // TODO [CL] Those files will be written to the working directory.
        //           Relative path needed.
        M.write("global_matrix_M.txt");
        K.write("global_matrix_K.txt");
        b.write("global_vector_b.txt");
    }
#endif

#if 0 && defined(OGS_USE_EIGENLIS)
    // TODO put that somewhere

    // scale back
    for (unsigned dof = 0; dof < NODAL_DOF; ++dof)
    {
        auto const& trafo = (dof == 0) ? _assembly_params.trafo_p
                          : (dof == 1) ? _assembly_params.trafo_T
                                       : _assembly_params.trafo_x;

        for (std::size_t i = 0; i < BP::_mesh.getNNodes(); ++i)
        {
            MeshLib::Location loc(BP::_mesh.getID(), MeshLib::MeshItemType::Node, i);
            auto const idx = BP::_local_to_global_index_map->getGlobalIndex(loc, dof);

            // TODO: _A
            x_curr[idx] *= trafo.dxdy(0);
        }
    }
#endif

    if (BP::_process_output.output_iteration_results)
    {
        DBUG("output results of iteration %li", _total_iteration);
        std::string fn = "tes_iter_" + std::to_string(_total_iteration) +
                         + "_ts_" + std::to_string(_timestep)
                         + "_" +    std::to_string(_assembly_params.iteration_in_current_timestep)
                         + "_" +    std::to_string(num_try)
                         + ".vtu";

        output(fn, 0);
    }

    ++num_try;
}

#endif
