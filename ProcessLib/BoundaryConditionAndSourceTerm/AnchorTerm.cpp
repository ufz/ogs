// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AnchorTerm.h"

#include <Eigen/Core>
#include <cassert>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"

namespace ProcessLib
{
template <int GlobalDim>
AnchorTerm<GlobalDim>::AnchorTerm(
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::size_t const source_term_mesh_id,
    MeshLib::Mesh const& st_mesh,
    const int variable_id,
    ParameterLib::Parameter<double> const& parameter)
    : SourceTerm(std::move(source_term_dof_table)),
      source_term_mesh_id_(source_term_mesh_id),
      st_mesh_(st_mesh),
      variable_id_(variable_id),
      parameter_(parameter)
{
    DBUG("Create AnchorTerm.");
}

template <int GlobalDim>
void AnchorTerm<GlobalDim>::integrate(const double t, GlobalVector const& x,
                                      GlobalVector& b, GlobalMatrix* jac) const
{
    DBUG("Assemble AnchorTerm.");

    using GlobalDimVector = Eigen::Vector<double, GlobalDim>;
    using GlobalDimMatrix = Eigen::Matrix<double, GlobalDim, GlobalDim>;

    for (MeshLib::Element const* const element : st_mesh_.getElements())
    {
        auto const element_id = element->getID();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(element_id);

        std::vector<GlobalIndexType> const global_indices =
            NumLib::getIndices(element_id, *_source_term_dof_table);
        assert(global_indices.size() == 2 * GlobalDim);

        Eigen::Vector<double, 2 * GlobalDim> const local_x =
            MathLib::toVector(x.get(global_indices));

        Eigen::Vector<double, 2 * GlobalDim> local_rhs =
            Eigen::Vector<double, 2 * GlobalDim>::Zero();
        Eigen::Matrix<double, 2 * GlobalDim, 2 * GlobalDim> local_Jac =
            Eigen::Matrix<double, 2 * GlobalDim, 2 * GlobalDim>::Zero();

        // The local indices are, due to the nature of the DOF table, all even
        // for the first node and odd for the second node.
        auto node_local_indices = [](int const i)
        { return Eigen::seqN(i, Eigen::fix<GlobalDim>, Eigen::fix<2>); };

        auto node_coords = [element](int const i)
        { return element->getNode(i)->asEigenVector3d(); };
        GlobalDimVector const l_original =
            (node_coords(1) - node_coords(0)).template head<GlobalDim>();
        double const l_original_norm = l_original.norm();

        // Displacement in the two nodes.
        auto u = [&local_x, &node_local_indices](int const i)
        { return local_x(node_local_indices(i)); };
        GlobalDimVector const l = l_original + u(1) - u(0);

        double const K = parameter_(t, pos)[0];
        GlobalDimVector const f = l_original / l_original_norm * K *
                                  (l.norm() - l_original_norm) /
                                  l_original_norm;

        GlobalDimMatrix const Df = l_original / l_original_norm * K *
                                   l.transpose() / l.norm() / l_original_norm;

        // Signs for the two nodes alternate.
        constexpr auto even_odd_sign = [](int const n)
        { return (n % 2 == 0) ? 1.0 : -1.0; };

        // There are two blocks in the rhs vector for the two nodes of the
        // element and four blocks in the local Jacobian formed by the four
        // combinations of the two nodes.
        for (int i = 0; i < 2; ++i)
        {
            local_rhs(node_local_indices(i)).noalias() += even_odd_sign(i) * f;

            for (int j = 0; j < 2; ++j)
            {
                local_Jac(node_local_indices(i), node_local_indices(j))
                    .noalias() += even_odd_sign(i) * even_odd_sign(j) * Df;
            }
        }

        b.add(global_indices, local_rhs);
        if (jac)
        {
            jac->add({global_indices, global_indices}, local_Jac);
        }
    }
}

template class AnchorTerm<2>;
template class AnchorTerm<3>;

}  // namespace ProcessLib
