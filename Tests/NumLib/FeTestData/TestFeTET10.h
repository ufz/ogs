/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 13, 2017, 2:57 PM
 */

#pragma once

#include "MeshLib/Elements/Tet.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
class TestFeTET10
{
public:
    using ShapeFunction = NumLib::ShapeTet10;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeTET10<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Tet10;
    static const unsigned dim = 3;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 5;
    static const unsigned n_sample_pt_order3 = 15;
    static const unsigned global_dim = MeshElementType::dimension;
    const double perturbation = 0.2;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        // nodes on the base
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.0, 1.0, 0.0);
        // top node
        nodes[3] = new MeshLib::Node(0.0, 0.0, 1.0);

        // mid-edge nodes on the base
        nodes[4] = new MeshLib::Node(0.5, 0.0, 0.0);
        nodes[5] = new MeshLib::Node(0.5, 0.5, 0.0);
        nodes[6] = new MeshLib::Node(0.0, 0.5, 0.0);

        // mid-edge nodes between the base and the top
        nodes[7] = new MeshLib::Node(perturbation, 0.0, 0.5);
        nodes[8] = new MeshLib::Node(0.5 + perturbation, 0.0, 0.5);
        nodes[9] = new MeshLib::Node(perturbation, 0.5, 0.5);

        return new MeshElementType(nodes);
    }

    double getVolume() const { return 1.0 / 6.; }
};

}  // namespace FeTestData
