/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Edge.h"

namespace MeshLib {

#ifndef WIN32
/// \todo Windows compiler does not accept this definition and issues a linking error.
const unsigned Edge::dimension;
#endif

} /* namespace */

