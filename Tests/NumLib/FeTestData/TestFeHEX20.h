/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 13, 2017, 2:40 PM
 */

#pragma once
#include "MeshLib/Elements/Hex.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
class TestFeHEX20
{
public:
    using ShapeFunction = NumLib::ShapeHex20;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeHEX20<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Hex20;
    static const unsigned dim = 3;  // MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 2 * 2 * 2;
    static const unsigned n_sample_pt_order3 = 3 * 3 * 3;
    static const unsigned global_dim = MeshElementType::dimension;
    const double perturbation = 0.2;

    /// create a 20 node hexahedral element with curved edges
    MeshElementType* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        // bottom
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(1.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(1.0, 1.0, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 1.0, 0.0);
        // top
        nodes[4] = new MeshLib::Node(0.0, 0.0, 1.0);
        nodes[5] = new MeshLib::Node(1.0, 0.0, 1.0);
        nodes[6] = new MeshLib::Node(1.0, 1.0, 1.0);
        nodes[7] = new MeshLib::Node(0.0, 1.0, 1.0);

        // mid-edge nodes on bottom
        nodes[8] = new MeshLib::Node(0.5, perturbation, 0.0);
        nodes[9] = new MeshLib::Node(1.0 + perturbation, 0.5, 0.0);
        nodes[10] = new MeshLib::Node(0.5, 1.0 + perturbation, 0.0);
        nodes[11] = new MeshLib::Node(perturbation, 0.5, 0.0);

        // mid-edge nodes on top
        nodes[12] = new MeshLib::Node(0.5, perturbation, 1.0);
        nodes[13] = new MeshLib::Node(1.0 + perturbation, 0.5, 1.0);
        nodes[14] = new MeshLib::Node(0.5, 1.0 + perturbation, 1.0);
        nodes[15] = new MeshLib::Node(perturbation, 0.5, 1.0);

        // mid-edge nodes between the bottom and the top
        nodes[16] = new MeshLib::Node(perturbation, 0.0, 0.5);
        nodes[17] = new MeshLib::Node(1.0 + perturbation, 0.0, 0.5);
        nodes[18] = new MeshLib::Node(1.0 + perturbation, 1.0, 0.5);
        nodes[19] = new MeshLib::Node(perturbation, 1.0, 0.5);

        return new MeshElementType(nodes);
    }

    double getVolume() const { return 1.0; }
};

}  // namespace FeTestData
