/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestFePYRA13.h
 *  Created on July 13, 2017, 2:57 PM
 */

#pragma once

#include "MeshLib/Elements/Pyramid.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
class TestFePYRA13
{
public:
    using ShapeFunction = NumLib::ShapePyra13;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FePYRA13<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Pyramid13;
    static const unsigned dim = 3;  // MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 5;
    static const unsigned n_sample_pt_order3 = 13;
    static const unsigned global_dim = MeshElementType::dimension;
    const double perturbation = 0.2;

    /// create a mesh element
    MeshElementType* createMeshElement()
    {
        auto** nodes = new MeshLib::Node*[e_nnodes];
        // nodes on base
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.5);
        nodes[1] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.0, 0.5, 0.0);
        nodes[3] = new MeshLib::Node(0.0, 0.5, 0.5);
        // node at the top
        nodes[4] = new MeshLib::Node(0.25, 0.25, 0.25);

        // mid-edge nodes on base
        nodes[5] = new MeshLib::Node(0.0, perturbation, 0.25);
        nodes[6] = new MeshLib::Node(0.0, 0.25, 0.0);
        nodes[7] = new MeshLib::Node(0.0, 0.5 + perturbation, 0.25);
        nodes[8] = new MeshLib::Node(0.0, 0.25, 0.5);

        // mid-edge nodes between the base and the top
        nodes[9] = new MeshLib::Node(0.125, 0.125, 0.375 + perturbation);
        nodes[10] = new MeshLib::Node(0.125, 0.125, 0.125 + perturbation);
        nodes[11] = new MeshLib::Node(0.125, 0.375, 0.125 + perturbation);
        nodes[12] = new MeshLib::Node(0.125, 0.375, 0.375 + perturbation);

        return new MeshElementType(nodes);
    }

    double getVolume() const { return 0.5 * 0.5 * 0.25 / 3.0; }
};

}  // namespace
