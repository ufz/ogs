/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "NumLib/NumericsConfig.h"

#include "GlobalMatrixProviders.h"
#include "SimpleMatrixVectorProvider.h"


// Initializes the static members of the structs in the header file
// associated with this file.
#define INITIALIZE_GLOBAL_MATRIX_VECTOR_PROVIDER(MAT, VEC, VARNAME) \
    static std::unique_ptr<MathLib::SimpleMatrixVectorProvider<MAT, VEC>> VARNAME{ \
        new MathLib::SimpleMatrixVectorProvider<MAT, VEC>}; \
    \
    namespace MathLib { \
    template<> \
    VectorProvider<VEC>& GlobalVectorProvider<VEC>::provider = *VARNAME; \
    \
    template<> \
    MatrixProvider<MAT>& GlobalMatrixProvider<MAT>::provider = *VARNAME; \
    }


#ifdef OGS_USE_EIGEN

INITIALIZE_GLOBAL_MATRIX_VECTOR_PROVIDER(Eigen::MatrixXd, Eigen::VectorXd,
                                         eigenGlobalMatrixVectorProvider)

#endif


using GlobalMatrix = detail::GlobalMatrixType;
using GlobalVector = detail::GlobalVectorType;

INITIALIZE_GLOBAL_MATRIX_VECTOR_PROVIDER(GlobalMatrix, GlobalVector,
                                         globalSetupGlobalMatrixVectorProvider)


namespace MathLib
{
void cleanupGlobalMatrixProviders()
{
    eigenGlobalMatrixVectorProvider.reset();
    globalSetupGlobalMatrixVectorProvider.reset();
}
}
