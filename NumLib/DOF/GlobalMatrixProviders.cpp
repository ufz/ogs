// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
