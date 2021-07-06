/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstddef>

#include "MathLib/LinAlg/MatrixSpecifications.h"

namespace NumLib
{
/*! Manages storage for matrices.
 *
 * This the matrix-analog of VectorProvider. The same notes apply to this class.
 */
class MatrixProvider
{
public:
    //! Get an uninitialized matrix with the given \c id.
    virtual GlobalMatrix& getMatrix(std::size_t& id) = 0;

    //! Get a matrix according to the given specifications in the storage
    //! of the matrix with the given \c id.
    virtual GlobalMatrix& getMatrix(MathLib::MatrixSpecifications const& ms,
                                    std::size_t& id) = 0;

    //! Release the given matrix.
    //!
    //! \pre \c A must have been acquired before, i.e., you must not call this
    //! method twice in a row in the same \c A!
    virtual void releaseMatrix(GlobalMatrix const& A) = 0;

    virtual ~MatrixProvider() = default;
};
}  // namespace NumLib
