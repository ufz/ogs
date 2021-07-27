/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
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
          _data(data),
          _local_matrix_size(local_matrix_size),
          _surface_normal(getOrientedSurfaceNormal(e))
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, std::vector<GlobalVector*> const& x,
                  int const process_id, GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType _local_rhs = NodalVectorType::Zero(_local_matrix_size);
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const boundary_permeability_node_values =
            _data.boundary_permeability.getNodalValuesOnElement(Base::_element,
                                                                t);
        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();

        auto const indices =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);
        std::vector<double> const local_values = x[process_id]->get(indices);
        std::size_t const bulk_element_id =
            _data.bulk_element_ids[Base::_element.getID()];
        std::size_t const bulk_face_id =
            _data.bulk_face_ids[Base::_element.getID()];
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;
            auto const& wp = Base::_integration_method.getWeightedPoint(ip);

            auto const bulk_element_point = MeshLib::getBulkElementPoint(
                _data.process.getMesh(), bulk_element_id, bulk_face_id, wp);

            double int_pt_value = 0.0;
            NumLib::shapeFunctionInterpolate(local_values, N, int_pt_value);

            NodalVectorType const neumann_node_values =
                -boundary_permeability_node_values * int_pt_value *
                _data.process.getFlux(bulk_element_id, bulk_element_point, t, x)
                    .dot(_surface_normal);
            _local_rhs.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices, _local_rhs);
    }

private:
    Eigen::Vector3d getOrientedSurfaceNormal(MeshLib::Element const& e) const
    {
        // At the moment (2016-09-28) the surface normal is not oriented
        // according to the right hand rule
        // for correct results it is necessary to multiply the normal with -1
        Eigen::Vector3d surface_normal =
            -MeshLib::FaceRule::getSurfaceNormal(e).normalized();
        auto const zeros_size = 3 - _data.process.getMesh().getDimension();
        surface_normal.tail(zeros_size).setZero();
        return surface_normal;
    }

    HCNonAdvectiveFreeComponentFlowBoundaryConditionData const& _data;
    std::size_t const _local_matrix_size;
    Eigen::Vector3d const _surface_normal;
};

}  // namespace ProcessLib
