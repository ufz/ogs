/**
 * \file ShapeStaticConsts.cpp
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The purpose of this file is to prevent linker errors.
 *
 * The classes defined in the included headers each have <tt>static const</tt> members
 * \c DIM and \c NPOINTS. The compiler is allowed to optimize those symbols out.
 * Therefore the linker later on might not find those symbols anymore and raises errors.
 * The <tt>static const</tt> members are explicitly instantiated here in order to rule
 * such errors out.
 *
 * Helper to generate that code:
 *
 * \code{.unparsed}
 * find . \( -name '*.h' -and -not -name '*-impl.h' \) -printf '%f\n' \
 * | sort \
 * | while read f; do
 *     c="${f%.h}";
 *     echo "const std::size_t $c::DIM;";
 *     echo "const std::size_t $c::NPOINTS;";
 * done >ShapeStaticConsts.cpp
 * \endcode
 *
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

const std::size_t ShapeHex20::DIM;
const std::size_t ShapeHex20::NPOINTS;
const std::size_t ShapeHex8::DIM;
const std::size_t ShapeHex8::NPOINTS;

const std::size_t ShapeLine2::DIM;
const std::size_t ShapeLine2::NPOINTS;
const std::size_t ShapeLine3::DIM;
const std::size_t ShapeLine3::NPOINTS;

const std::size_t ShapePoint1::DIM;
const std::size_t ShapePoint1::NPOINTS;

const std::size_t ShapePrism15::DIM;
const std::size_t ShapePrism15::NPOINTS;
const std::size_t ShapePrism6::DIM;
const std::size_t ShapePrism6::NPOINTS;

const std::size_t ShapePyra13::DIM;
const std::size_t ShapePyra13::NPOINTS;
const std::size_t ShapePyra5::DIM;
const std::size_t ShapePyra5::NPOINTS;

const std::size_t ShapeQuad4::DIM;
const std::size_t ShapeQuad4::NPOINTS;
const std::size_t ShapeQuad8::DIM;
const std::size_t ShapeQuad8::NPOINTS;
const std::size_t ShapeQuad9::DIM;
const std::size_t ShapeQuad9::NPOINTS;

const std::size_t ShapeTet10::DIM;
const std::size_t ShapeTet10::NPOINTS;
const std::size_t ShapeTet4::DIM;
const std::size_t ShapeTet4::NPOINTS;

const std::size_t ShapeTri3::DIM;
const std::size_t ShapeTri3::NPOINTS;
const std::size_t ShapeTri6::DIM;
const std::size_t ShapeTri6::NPOINTS;

}
