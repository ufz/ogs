/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MatrixProviderUser.h"
#include "MatrixVectorTraits.h"

#ifdef OGS_USE_EIGEN

namespace MathLib
{

std::unique_ptr<Eigen::MatrixXd>
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance()
{
    return std::unique_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd);
}

std::unique_ptr<Eigen::MatrixXd>
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance(Eigen::MatrixXd const& A)
{
    return std::unique_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(A));
}

std::unique_ptr<Eigen::MatrixXd>
MatrixVectorTraits<Eigen::MatrixXd>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows = spec.dof_table ? spec.dof_table->dofSizeWithGhosts() : spec.nrows;
    auto const ncols = spec.dof_table ? nrows : spec.ncols;

    return std::unique_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(nrows, ncols));
}

std::unique_ptr<Eigen::VectorXd>
MatrixVectorTraits<Eigen::VectorXd>::
newInstance()
{
    return std::unique_ptr<Eigen::VectorXd>(new Eigen::VectorXd);
}

std::unique_ptr<Eigen::VectorXd>
MatrixVectorTraits<Eigen::VectorXd>::
newInstance(Eigen::VectorXd const& A)
{
    return std::unique_ptr<Eigen::VectorXd>(new Eigen::VectorXd(A));
}

std::unique_ptr<Eigen::VectorXd>
MatrixVectorTraits<Eigen::VectorXd>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows = spec.dof_table ? spec.dof_table->dofSizeWithGhosts() : spec.nrows;

    return std::unique_ptr<Eigen::VectorXd>(new Eigen::VectorXd(nrows));
}

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

#include "MeshLib/NodePartitionedMesh.h"

namespace MathLib
{

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance()
{
    return std::unique_ptr<PETScMatrix>(new PETScMatrix);
}

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance(PETScMatrix const& A)
{
    return std::unique_ptr<PETScMatrix>(new PETScMatrix(A));
}

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows =
        spec.dof_table ? spec.dof_table->dofSizeWithoutGhosts() : spec.nrows;
    auto const ncols = spec.dof_table ? nrows : spec.ncols;

    // TODO I guess it is not hard to make AssemblerLib::computeSparsityPattern()
    //      work also for NodePartitionedMesh'es.
    assert(((!spec.sparsity_pattern) || spec.sparsity_pattern->size() == 0u) &&
           "The sparsity pattern is not used with PETSc so I rather crash than"
           " silently accept a precomputed sparsity pattern.");

    if (spec.mesh)
    {
        assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(spec.mesh));
        auto const& mesh = *static_cast<MeshLib::NodePartitionedMesh const*>(spec.mesh);

        auto const max_nonzeroes = mesh.getMaximumNConnectedNodesToNode();

        PETScMatrixOption mat_opt;
        mat_opt.d_nz = max_nonzeroes;
        mat_opt.o_nz = max_nonzeroes;
        mat_opt.is_global_size = false;
        return std::unique_ptr<PETScMatrix>(
            new PETScMatrix(nrows, ncols, mat_opt));
    }
    else
        return std::unique_ptr<PETScMatrix>(
            new PETScMatrix(nrows, ncols));
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance()
{
    return std::unique_ptr<PETScVector>(new PETScVector);
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance(PETScVector const& x)
{
    return std::unique_ptr<PETScVector>(new PETScVector(x));
}

std::unique_ptr<PETScVector>
MatrixVectorTraits<PETScVector>::
newInstance(MatrixSpecifications const& spec)
{
    auto const is_global_size = false;

    if (spec.dof_table) {
        auto const& dt = *spec.dof_table;
        return std::unique_ptr<PETScVector>(new PETScVector(
            dt.dofSizeWithoutGhosts(), dt.getGhostIndices(), is_global_size));
    } else {
        return std::unique_ptr<PETScVector>(
            new PETScVector(spec.nrows, is_global_size));
    }
}

} // namespace MathLib


#elif defined(OGS_USE_EIGEN)

namespace MathLib
{

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance()
{
    return std::unique_ptr<EigenMatrix>(new EigenMatrix(0, 0)); // TODO default constructor
}

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance(EigenMatrix const& A)
{
    return std::unique_ptr<EigenMatrix>(new EigenMatrix(A));
}

std::unique_ptr<EigenMatrix>
MatrixVectorTraits<EigenMatrix>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows = spec.dof_table ? spec.dof_table->dofSizeWithGhosts() : spec.nrows;

    auto A = std::unique_ptr<EigenMatrix>(new EigenMatrix(nrows));

    if (spec.sparsity_pattern)
        setMatrixSparsity(*A, *spec.sparsity_pattern);

    return A;
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance()
{
    return std::unique_ptr<EigenVector>(new EigenVector);
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance(EigenVector const& x)
{
    return std::unique_ptr<EigenVector>(new EigenVector(x));
}

std::unique_ptr<EigenVector>
MatrixVectorTraits<EigenVector>::
newInstance(MatrixSpecifications const& spec)
{
    auto const nrows = spec.dof_table ? spec.dof_table->dofSizeWithGhosts() : spec.nrows;

    return std::unique_ptr<EigenVector>(new EigenVector(nrows));
}

} // namespace MathLib

#endif // defined(OGS_USE_EIGEN)
