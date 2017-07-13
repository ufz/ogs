/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestFeQuad9.h
 *  Created on July 13, 2017, 2:52 PM
 */

#pragma once
#include "MeshLib/Elements/Quad.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
class TestFeQUAD9
{
public:
    using ShapeFunction = NumLib::ShapeQuad9;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeQUAD9<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Quad9;
    static const unsigned dim = 2;  // MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2 * 2;
    static const unsigned n_sample_pt_order3 = 3 * 3;
    static const unsigned global_dim = MeshElementType::dimension;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        // square
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);

        nodes[4] = new MeshLib::Node(0.5, 0.0, 0.0);
        nodes[5] = new MeshLib::Node(1.0, 0.5, 0.0);
        nodes[6] = new MeshLib::Node(0.5, 1.0, 0.0);
        nodes[7] = new MeshLib::Node(0.0, 0.5, 0.0);

        nodes[8] = new MeshLib::Node(0.5, 0.5, 0.0);
        return new MeshElementType(nodes);
    }

    double getVolume() const { return 1.0; }
};

}  // namespace
