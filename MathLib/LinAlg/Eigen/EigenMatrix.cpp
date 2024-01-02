/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
