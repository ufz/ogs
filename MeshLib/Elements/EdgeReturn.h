// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

namespace MeshLib
{

class Element;

/// Returns always null pointer
class NoEdgeReturn
{
public:
    /// Returns i-th edge of the given element
    static const Element* getEdge(const Element* /*e*/, unsigned /*i*/)
    {
        return nullptr;
    }
};

/// Returns linear order edge
class LinearEdgeReturn
{
public:
    /// Returns i-th edge of the given element
    static const Element* getEdge(const Element* e, unsigned i);
};

/// Returns quadratic order edge
class QuadraticEdgeReturn
{
public:
    /// Returns i-th edge of the given element
    static const Element* getEdge(const Element* e, unsigned i);
};

}  // namespace MeshLib
