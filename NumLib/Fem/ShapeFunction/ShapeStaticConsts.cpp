/**
 * \file ShapeStaticConsts.cpp
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The purpose of this file is to prevent linker errors.
 *
 * The classes defined in the included headers each have <tt>static const</tt>
 * members \c DIM and \c NPOINTS. The compiler is allowed to optimize those
 * symbols out. Therefore the linker later on might not find those symbols
 * anymore and raises errors. The <tt>static const</tt> members are explicitly
 * instantiated here in order to rule such errors out.
 *
 * Helper to generate that code:
 *
 * \code{.sh}
 * find . '(' -name '*.h' -and -not -name '*-impl.h' ')' -printf '%f\n' \
 * | sort \
 * | while read f; do
 *     c="${f%.h}";
 *     echo "const unsigned $c::DIM;";
 *     echo "const unsigned $c::NPOINTS;";
 * done >ShapeStaticConsts.cpp
 * \endcode
 */

#include "ShapeHex20.h"
#include "ShapeHex8.h"
#include "ShapeLine2.h"
#include "ShapeLine3.h"
#include "ShapePoint1.h"
#include "ShapePrism15.h"
#include "ShapePrism6.h"
#include "ShapePyra13.h"
#include "ShapePyra5.h"
#include "ShapeQuad4.h"
#include "ShapeQuad8.h"
#include "ShapeQuad9.h"
#include "ShapeTet10.h"
#include "ShapeTet4.h"
#include "ShapeTri3.h"
#include "ShapeTri6.h"

namespace NumLib
{

const unsigned ShapeHex20::DIM;
const unsigned ShapeHex20::NPOINTS;
const unsigned ShapeHex8::DIM;
const unsigned ShapeHex8::NPOINTS;

const unsigned ShapeLine2::DIM;
const unsigned ShapeLine2::NPOINTS;
const unsigned ShapeLine3::DIM;
const unsigned ShapeLine3::NPOINTS;

const unsigned ShapePoint1::DIM;
const unsigned ShapePoint1::NPOINTS;

const unsigned ShapePrism15::DIM;
const unsigned ShapePrism15::NPOINTS;
const unsigned ShapePrism6::DIM;
const unsigned ShapePrism6::NPOINTS;

const unsigned ShapePyra13::DIM;
const unsigned ShapePyra13::NPOINTS;
const unsigned ShapePyra5::DIM;
const unsigned ShapePyra5::NPOINTS;

const unsigned ShapeQuad4::DIM;
const unsigned ShapeQuad4::NPOINTS;
const unsigned ShapeQuad8::DIM;
const unsigned ShapeQuad8::NPOINTS;
const unsigned ShapeQuad9::DIM;
const unsigned ShapeQuad9::NPOINTS;

const unsigned ShapeTet10::DIM;
const unsigned ShapeTet10::NPOINTS;
const unsigned ShapeTet4::DIM;
const unsigned ShapeTet4::NPOINTS;

const unsigned ShapeTri3::DIM;
const unsigned ShapeTri3::NPOINTS;
const unsigned ShapeTri6::DIM;
const unsigned ShapeTri6::NPOINTS;

}
