/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GlobalMatrixProviders.h"

#include <memory>

#include "SimpleMatrixVectorProvider.h"

static std::unique_ptr<NumLib::SimpleMatrixVectorProvider>
    globalSetupGlobalMatrixVectorProvider =
        std::make_unique<NumLib::SimpleMatrixVectorProvider>();

namespace NumLib
{
VectorProvider& GlobalVectorProvider::provider =
    *(globalSetupGlobalMatrixVectorProvider);

MatrixProvider& GlobalMatrixProvider::provider =
    *(globalSetupGlobalMatrixVectorProvider);

void cleanupGlobalMatrixProviders()
{
    globalSetupGlobalMatrixVectorProvider->clear();
}
}  // namespace NumLib
