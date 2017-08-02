/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestFePRISM15.h
 *  Created on July 13, 2017, 2:56 PM
 */

#pragma once
#include "MeshLib/Elements/Prism.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
class TestFePRISM15
{
public:
    using ShapeFunction = NumLib::ShapePrism15;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FePRISM15<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Prism15;
    static const unsigned dim = 3;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 6;
    static const unsigned n_sample_pt_order3 = 6;  // TODO no implementation yet
    static const unsigned global_dim = MeshElementType::dimension;
    const double perturbation = 0.2;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        // vertex nodes on bottom
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(0.5, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.5, 0.5, 0.0);

        // vertex nodes on top
        nodes[3] = new MeshLib::Node(0.0, 0.0, 0.5);
        nodes[4] = new MeshLib::Node(0.5, 0.0, 0.5);
        nodes[5] = new MeshLib::Node(0.5, 0.5, 0.5);

        // mid-edge nodes on bottom
        nodes[6] = new MeshLib::Node(0.25, 0.0, perturbation);
        nodes[7] = new MeshLib::Node(0.5, 0.25, perturbation);
        nodes[8] = new MeshLib::Node(0.25, 0.25, perturbation);

        // mid-edge nodes on top
        nodes[9] = new MeshLib::Node(0.25, 0.0, 0.5 + perturbation);
        nodes[10] = new MeshLib::Node(0.5, 0.25, 0.5 + perturbation);
        nodes[11] = new MeshLib::Node(0.25, 0.25, 0.5 + perturbation);

        // mid-edge nodes between the bottom and the top
        nodes[12] = new MeshLib::Node(0.0, perturbation, 0.25);
        nodes[13] = new MeshLib::Node(0.5, perturbation, 0.25);
        nodes[14] = new MeshLib::Node(0.5, 0.5 + perturbation, 0.25);

        return new MeshElementType(nodes);
    }

    double getVolume() const { return 0.0625; }
};

}  // namespace
