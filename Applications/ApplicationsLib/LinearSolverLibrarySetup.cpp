// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "LinearSolverLibrarySetup.h"

#include <memory>

#include "BaseLib/Logging.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"

// All concrete LinearSolverLibrarySetup implementations are in the detail
// namespace. They should not be instantiated manually, but only via
// LinearSolverLibrarySetup::create().
#if defined(USE_PETSC)
#include <mpi.h>
#include <petsc.h>

#include "BaseLib/MPI.h"

namespace ApplicationsLib::detail
{
struct LinearSolverLibrarySetupPETSc final
    : public ApplicationsLib::LinearSolverLibrarySetup
{
    LinearSolverLibrarySetupPETSc(int argc, char* argv[])
    {
        char help[] = "ogs6 with PETSc \n";
        PETSC_COMM_WORLD = BaseLib::MPI::OGS_COMM_WORLD;
        PetscInitialize(&argc, &argv, nullptr, help);
        MPI_Comm_set_errhandler(PETSC_COMM_WORLD, MPI_ERRORS_RETURN);
    }

    ~LinearSolverLibrarySetupPETSc()
    {
        NumLib::cleanupGlobalMatrixProviders();
        PetscFinalize();
    }
};
using LinearSolverLibrarySetupImpl = LinearSolverLibrarySetupPETSc;
}  // namespace ApplicationsLib::detail
#elif defined(USE_LIS)
#include "MathLib/LinAlg/Lis/LisWrapper.h"
namespace ApplicationsLib::detail
{
struct LinearSolverLibrarySetupLis final
    : public ApplicationsLib::LinearSolverLibrarySetup
{
    LinearSolverLibrarySetupLis(int argc, char* argv[])
    {
        lis_initialize(&argc, &argv);
    }

    ~LinearSolverLibrarySetupLis()
    {
        NumLib::cleanupGlobalMatrixProviders();
        lis_finalize();
    }
};
using LinearSolverLibrarySetupImpl = LinearSolverLibrarySetupLis;
}  // namespace ApplicationsLib::detail
#else
namespace ApplicationsLib::detail
{
struct LinearSolverLibrarySetupEigen final
    : public ApplicationsLib::LinearSolverLibrarySetup
{
    LinearSolverLibrarySetupEigen(int /*argc*/, char* /*argv*/[])
    {
#ifdef _OPENMP
        const char* omp_num_threads_env = std::getenv("OMP_NUM_THREADS");
        if (omp_num_threads_env)
        {
            INFO("OMP_NUM_THREADS is set to: {:s}", omp_num_threads_env);
        }
        else
        {
            WARN("OMP_NUM_THREADS is not set, falling back to: {:d}",
                 omp_get_max_threads());
        }
#endif
        INFO("Eigen use {:d} threads", Eigen::nbThreads());
    }
    ~LinearSolverLibrarySetupEigen() override
    {
        NumLib::cleanupGlobalMatrixProviders();
    }
};
using LinearSolverLibrarySetupImpl = LinearSolverLibrarySetupEigen;
}  // namespace ApplicationsLib::detail
#endif

namespace ApplicationsLib
{
LinearSolverLibrarySetup::~LinearSolverLibrarySetup()
{
    DBUG("Tearing down linear solver library setup.");
}

std::shared_ptr<LinearSolverLibrarySetup> LinearSolverLibrarySetup::create(
    int argc, char* argv[])
{
    static std::weak_ptr<LinearSolverLibrarySetup> instance;
    static std::mutex mutex;

    // Lock to avoid multiple initializations.
    std::lock_guard lock{mutex};

    std::shared_ptr lsls = instance.lock();

    if (!lsls)
    {
        DBUG("Initializing linear solver library for the first time");
        lsls =
            std::make_shared<detail::LinearSolverLibrarySetupImpl>(argc, argv);

        // Necessary such that the internally cached ptr and the returned ptr
        // both point to the same object, which might be returned from
        // subsequent create() calls.
        instance = lsls;
    }
    else
    {
        DBUG("Initializing linear solver library for the {}th time",
             lsls.use_count());
#if defined(USE_PETSC) || defined(USE_LIS)
        // PETSc and Lis parse commandline arguments. Therefore, we exclude them
        // for the time being such that nobody can accidentally pass
        // different arguments to different linear solver library setups.
        OGS_FATAL(
            "Reusing a LinearSolverLibrarySetup has not been defined for "
            "PETSc or LIS. In the present build configuration you cannot run "
            "two separate OGS simulations in the same process.");
#else
        // No initialization needed for Eigen, just return the internally cached
        // ptr.
#endif
    }

    // Consistency check that the internally cached ptr and the returned ptr
    // both point to the same object.
    // C++26 will have owner_equal()
    assert(!instance.owner_before(lsls) && !lsls.owner_before(instance));

    return lsls;
}
}  // namespace ApplicationsLib
