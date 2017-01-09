/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_MATRIX_PROVIDER_USER_H
#define NUMLIB_MATRIX_PROVIDER_USER_H

#include <cstddef>

#include "MathLib/LinAlg/MatrixSpecifications.h"

namespace NumLib
{

class MatrixSpecificationsProvider
{
public:
    virtual MathLib::MatrixSpecifications getMatrixSpecifications() const = 0;

    virtual ~MatrixSpecificationsProvider() = default;
};


/*! Manages storage for vectors.
 *
 * This interface provides storage management semantics for vectors, which
 * can be acquired at a certain point in time, be released later on and
 * be acquired again etc.
 *
 * "Released" means that the caller indicates, that he currently does not need
 * that \c Vector instance, i.e, Either it is not needed anymore for the entire
 * program run, or it is temporarily not needed, but might be aquired again
 * later on. Thereby the implementation of this interface can decide, whether
 * the storage for the specific vector can be freed or reused in the meantime.
 *
 * All get-methods of this class come in two variants: One with an \c id argument,
 * and one without. The latter makes it possible to temporarily release a vector
 * during the program run and re-acquire it again. The former variant is intended
 * as a short-cut that simplifies vector acquisition for callers that will use
 * the same vector throughout all of their lifetime without releasing it
 * intermediately.
 *
 * \attention
 * The first time a vector is acquired by a method with \c id argument, the \c id
 * has to be initialized with zero. The respective method then will set the \c id
 * to the specific \c id of the returned vector.
 */
class VectorProvider
{
public:
    //! Get an uninitialized vector.
    virtual GlobalVector& getVector() = 0;

    //! Get an uninitialized vector with the given \c id.
    virtual GlobalVector& getVector(std::size_t& id) = 0;

    //! Get a copy of \c x.
    virtual GlobalVector& getVector(GlobalVector const& x) = 0;

    //! Get a copy of \c x in the storage of the vector with the given \c id.
    virtual GlobalVector& getVector(GlobalVector const& x, std::size_t& id) = 0;

    //! Get a vector according to the given specifications.
    virtual GlobalVector& getVector(MathLib::MatrixSpecifications const& ms) = 0;

    //! Get a vector according to the given specifications in the storage
    //! of the vector with the given \c id.
    virtual GlobalVector& getVector(MathLib::MatrixSpecifications const& ms, std::size_t& id) = 0;

    //! Release the given vector.
    //!
    //! \pre \c x must have been acquired before, i.e., you must not call this
    //! method twice in a row in the same \c x!
    virtual void releaseVector(GlobalVector const& x) = 0;

    virtual ~VectorProvider() = default;
};

/*! Manages storage for matrices.
 *
 * This the matrix-analog of VectorProvider. The same notes apply to this class.
 */
class MatrixProvider
{
public:
    //! Get an uninitialized matrix.
    virtual GlobalMatrix& getMatrix() = 0;

    //! Get an uninitialized matrix with the given \c id.
    virtual GlobalMatrix& getMatrix(std::size_t& id) = 0;

    //! Get a copy of \c A.
    virtual GlobalMatrix& getMatrix(GlobalMatrix const& A) = 0;

    //! Get a copy of \c A in the storage of the matrix with the given \c id.
    virtual GlobalMatrix& getMatrix(GlobalMatrix const& A, std::size_t& id) = 0;

    //! Get a matrix according to the given specifications.
    virtual GlobalMatrix& getMatrix(MathLib::MatrixSpecifications const& ms) = 0;

    //! Get a matrix according to the given specifications in the storage
    //! of the matrix with the given \c id.
    virtual GlobalMatrix& getMatrix(MathLib::MatrixSpecifications const& ms, std::size_t& id) = 0;

    //! Release the given matrix.
    //!
    //! \pre \c A must have been acquired before, i.e., you must not call this
    //! method twice in a row in the same \c A!
    virtual void releaseMatrix(GlobalMatrix const& A) = 0;

    virtual ~MatrixProvider() = default;
};

} // namespace NumLib

#endif // NUMLIB_MATRIX_PROVIDER_USER_H
