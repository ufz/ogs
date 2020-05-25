/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MeshLib/Elements/MapBulkElementPoint.h"
#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/MeshNodeParameter.h"
#include "ProcessLib/Process.h"

namespace ProcessLib
{
struct HCNonAdvectiveFreeComponentFlowBoundaryConditionData
{
    ParameterLib::Parameter<double> const& boundary_permeability;
    MeshLib::PropertyVector<std::size_t> const bulk_face_ids;
    MeshLib::PropertyVector<std::size_t> const bulk_element_ids;
    Process const& process;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;
    using NodalMatrixType = typename Base::NodalMatrixType;

public:
    /// The neumann_bc_term factor is directly integrated into the local
    /// element matrix.
    HCNonAdvectiveFreeComponentFlowBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HCNonAdvectiveFreeComponentFlowBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          data_(data),
          local_matrix_size_(local_matrix_size),
          surface_normal_(getOrientedSurfaceNormal(e))
    {
    }
    Eigen::RowVector3d getOrientedSurfaceNormal(MeshLib::Element const& e)
    {
        auto surface_element_normal = MeshLib::FaceRule::getSurfaceNormal(&e);
        surface_element_normal.normalize();
        // At the moment (2016-09-28) the surface normal is not oriented
        // according to the right hand rule
        // for correct results it is necessary to multiply the normal with -1
        surface_element_normal *= -1;
        return Eigen::Map<Eigen::RowVector3d const>(
            surface_element_normal.getCoords(),
            data_.process.getMesh().getDimension());
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& x,
                  int const process_id, GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType local_rhs_ = NodalVectorType::Zero(local_matrix_size_);
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const boundary_permeability_node_values =
            data_.boundary_permeability.getNodalValuesOnElement(Base::element_,
                                                                t);
        unsigned const n_integration_points =
            Base::integration_method_.getNumberOfPoints();

        auto const indices =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        std::vector<double> const local_values = x[process_id]->get(indices);
        std::size_t const bulk_element_id =
            data_.bulk_element_ids[Base::element_.getID()];
        std::size_t const bulk_face_id =
            data_.bulk_face_ids[Base::element_.getID()];
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::ns_and_weights_[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;
            auto const& wp = Base::integration_method_.getWeightedPoint(ip);

            auto const bulk_element_point = MeshLib::getBulkElementPoint(
                data_.process.getMesh(), bulk_element_id, bulk_face_id, wp);

            double int_pt_value = 0.0;
            NumLib::shapeFunctionInterpolate(local_values, N, int_pt_value);

            NodalVectorType const neumann_node_values =
                -boundary_permeability_node_values * int_pt_value *
                data_.process.getFlux(bulk_element_id, bulk_element_point, t, x)
                    .dot(surface_normal_);
            local_rhs_.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices, local_rhs_);
    }

private:
    HCNonAdvectiveFreeComponentFlowBoundaryConditionData const& data_;
    std::size_t const local_matrix_size_;
    Eigen::RowVector3d const surface_normal_;
};

}  // namespace ProcessLib
