// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EigenMatrix.h"

#include <fstream>

namespace MathLib
{

/// printout this matrix for debugging
void EigenMatrix::write(const std::string& filename) const
{
    std::ofstream of(filename);
    if (of)
    {
        write(of);
    }
}

/// printout this matrix for debugging
void EigenMatrix::write(std::ostream& os) const
{
    for (int k = 0; k < mat_.outerSize(); ++k)
    {
        for (RawMatrixType::InnerIterator it(mat_, k); it; ++it)
        {
            os << it.row() << " " << it.col() << " " << it.value() << "\n";
        }
    }
    os << std::endl;
}
}  // namespace MathLib
