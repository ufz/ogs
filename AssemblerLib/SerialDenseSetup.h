/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_SERIALDENSESETUP_H_
#define ASSEMBLERLIB_SERIALDENSESETUP_H_


#include "AssemblerLib/GlobalSetup.h"

#include "AssemblerLib/SerialDenseVectorMatrixBuilder.h"
#include "AssemblerLib/SerialExecutor.h"

namespace AssemblerLib
{

/// Using GlobalDenseMatrix and DenseVector for global entities and serial
/// global assembly loop.
typedef GlobalSetup<
        AssemblerLib::SerialDenseVectorMatrixBuilder,
        AssemblerLib::SerialExecutor>
    SerialDenseSetup;

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_SERIALDENSESETUP_H_
