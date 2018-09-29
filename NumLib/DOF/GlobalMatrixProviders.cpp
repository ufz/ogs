/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GlobalMatrixProviders.h"

#include <memory>

#include "SimpleMatrixVectorProvider.h"


// Initializes the static members of the structs in the header file
// associated with this file.
#define INITIALIZE_GLOBAL_MATRIX_VECTOR_PROVIDER(VARNAME) \
    static std::unique_ptr<NumLib::SimpleMatrixVectorProvider> VARNAME{ \
        new NumLib::SimpleMatrixVectorProvider}; \
    \
    namespace NumLib { \
    VectorProvider& GlobalVectorProvider::provider = *(VARNAME); \
    \
    MatrixProvider& GlobalMatrixProvider::provider = *(VARNAME); \
    }

INITIALIZE_GLOBAL_MATRIX_VECTOR_PROVIDER(globalSetupGlobalMatrixVectorProvider)


namespace NumLib
{
void cleanupGlobalMatrixProviders()
{
    globalSetupGlobalMatrixVectorProvider.reset();
}
}
