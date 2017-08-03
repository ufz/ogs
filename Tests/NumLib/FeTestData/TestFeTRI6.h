/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file   TestFeTRI6.h
 *  Created on July 13, 2017, 2:58 PM
 */

#pragma once

#include "MeshLib/Elements/Tri.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"

namespace FeTestData
{
/** Test the six node iso-parametric triangle element
 * Assuming that the edge points are translated outside of the square shape
 * element with a distance of b (perturbation) to their owner edges, we can
 * obtain the the area of the changed element as
 * \f[
 *   A = A_0 + \sum_{i=1}^3\frac{2}{3}a_i b
 * \f]
 * where\f$A_0\f$ is the area of the original triangle, \f$a_i\f$ is the length
 * of edge \f$i\f$ of the original element.
 */
class TestFeTRI6
{
public:
    using ShapeFunction = NumLib::ShapeTri6;

    // Fe type information
    template <template <typename> class ShapeMatrixPolicy_>
    using FeType = NumLib::FeTRI6<ShapeMatrixPolicy_>;

    using MeshElementType = MeshLib::Tri6;
    static const unsigned dim = 2;  // MeshElementType::dimension;
    static const unsigned e_nnodes = MeshElementType::n_all_nodes;
    static const unsigned n_sample_pt_order2 = 3;
    static const unsigned n_sample_pt_order3 = 4;
    static const unsigned global_dim = MeshElementType::dimension;

    const double a = 1.0;
    const double perturbation = 0.1;

    /// create a 6 node triangle element
    MeshElementType* createMeshElement()
    {
        // Concave 6 node triangle element with curved edges.
        auto** nodes = new MeshLib::Node*[e_nnodes];
        nodes[0] = new MeshLib::Node(0.0, 0.0, 0.0);
        nodes[1] = new MeshLib::Node(a, 0.0, 0.0);
        nodes[2] = new MeshLib::Node(0.0, a, 0.0);

        nodes[3] = new MeshLib::Node(0.5 * a, perturbation, 0.0);
        const double perturbation_hypotenuse =
            0.5 * perturbation * std::sqrt(2 * a * a);
        nodes[4] = new MeshLib::Node(0.5 * a - perturbation_hypotenuse,
                                     0.5 * a - perturbation_hypotenuse, 0.0);
        nodes[5] = new MeshLib::Node(perturbation, 0.5 * a, 0.0);

        return new MeshElementType(nodes);
    }

    double getVolume() const
    {
        // The length of hypotenuse is a_h sqrt(2 * a * a)
        // Area = a^2 /2 - 2 * 2 * a * b /3 - 2 * a_h * b /3
        // where b=perturbation
        return 0.5 * a * a - 4.0 * a * perturbation / 3.0 -
               2.0 * std::sqrt(2 * a * a) * perturbation / 3.0;
    }
};

}  // namespace
