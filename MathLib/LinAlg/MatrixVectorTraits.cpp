/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    return std::unique_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(spec.nrows, spec.ncols));
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
    return std::unique_ptr<Eigen::VectorXd>(new Eigen::VectorXd(spec.nrows));
}

} // namespace MathLib

#endif // OGS_USE_EIGEN


#ifdef USE_PETSC

namespace MathLib
{

std::unique_ptr<PETScMatrix>
MatrixVectorTraits<PETScMatrix>::
newInstance()
{
    return std::unique_ptr<PETScMatrix>(new PETScMatrix(0, 0)); // TODO default constructor
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
    return std::unique_ptr<PETScMatrix>(new PETScMatrix(spec.nrows)); // TODO sparsity pattern
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
    return std::unique_ptr<PETScVector>(new PETScVector(spec.nrows));
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
    return std::unique_ptr<EigenMatrix>(new EigenMatrix(spec.nrows)); // TODO sparsity pattern
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
    return std::unique_ptr<EigenVector>(new EigenVector(spec.nrows));
}

} // namespace MathLib

#endif // defined(OGS_USE_EIGEN)
