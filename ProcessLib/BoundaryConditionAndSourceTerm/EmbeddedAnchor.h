/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SourceTerm.h"

namespace ProcessLib
{
template <int GlobalDim>
class EmbeddedAnchor final : public SourceTermBase
{
public:
    explicit EmbeddedAnchor(MeshLib::Mesh const& bulk_mesh,
                            NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
                            std::size_t const source_term_mesh_id,
                            MeshLib::Mesh const& st_mesh, const int variable_id,
                            ParameterLib::Parameter<double> const& parameter);

    void getShapeMatricesAndGlobalIndicesAndDisplacements(
        MeshLib::Element const* const anchor_element,
        std::array<std::size_t, 2>& nodes_per_element,
        std::vector<Eigen::RowVectorXd>& shape_matrices,
        std::vector<GlobalIndexType>& global_indices,
        Eigen::Vector<double, 2 * GlobalDim>& local_x,
        GlobalVector const& x,
        ParameterLib::SpatialPosition& pos) const;

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk_;
    MeshLib::Mesh const& bulk_mesh_;
    std::size_t const source_term_mesh_id_;
    MeshLib::Mesh const& st_mesh_;
    int const variable_id_;
    std::array<int, GlobalDim> const component_ids_;
    ParameterLib::Parameter<double> const& parameter_;
    MeshLib::PropertyVector<std::size_t> const* bulk_element_ids_ = nullptr;
    MeshLib::PropertyVector<double> const* natural_coordinates_ = nullptr;
};

extern template class EmbeddedAnchor<2>;
extern template class EmbeddedAnchor<3>;
}  // namespace ProcessLib
