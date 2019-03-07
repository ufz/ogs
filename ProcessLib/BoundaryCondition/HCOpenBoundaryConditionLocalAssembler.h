/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Parameter/MeshNodeParameter.h"
#include "ProcessLib/Process.h"
#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "MeshLib/Elements/MapBulkElementPoint.h"

namespace ProcessLib
{
struct HCOpenBoundaryConditionData
{
    Parameter<double> const& boundary_permeability;
    MeshLib::PropertyVector<std::size_t> const bulk_face_ids;
    MeshLib::PropertyVector<std::size_t> const bulk_element_ids;
    Process const& process;

};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class HCOpenBoundaryConditionLocalAssembler final
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
    HCOpenBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        HCOpenBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_matrix_size(local_matrix_size)
    {
    }

    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* /*Jac*/) override
    {
        NodalVectorType _local_rhs(_local_matrix_size);
        _local_rhs.setZero();
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType const boundary_permeability_node_values =
            _data.boundary_permeability.getNodalValuesOnElement(Base::_element, t);
        auto surface_element_normal =
            MeshLib::FaceRule::getSurfaceNormal(&(Base::_element));
        surface_element_normal.normalize();
        // At the moment (2016-09-28) the surface normal is not oriented
        // according to the right hand rule
        // for correct results it is necessary to multiply the normal with -1
        surface_element_normal *= -1;
        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();
        auto const indices =
            NumLib::getIndices(mesh_item_id, dof_table_boundary);

        std::vector<double> const local_values =
            x.get(indices);

        std::size_t bulk_element_id = _data.bulk_element_ids[Base::_element.getID()];
        std::size_t bulk_face_id = _data.bulk_face_ids[Base::_element.getID()];
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        { 
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;

            auto const& wp = Base::_integration_method.getWeightedPoint(ip);

            auto const bulk_element_point = MeshLib::getBulkElementPoint(
                _data.process.getMesh(), bulk_element_id, bulk_face_id, wp);

            double int_pt_value = 0.0;

            NumLib::shapeFunctionInterpolate(local_values, N,
                                             int_pt_value);

            NodalVectorType const neumann_node_values = -boundary_permeability_node_values * int_pt_value *_data.process.getFlux(bulk_element_id, bulk_element_point, t, x).dot(Eigen::Map<Eigen::RowVectorXd const>(
                            surface_element_normal.getCoords(), _data.process.getMesh().getDimension()));
            _local_rhs.noalias() += N * neumann_node_values.dot(N) * w;
        }

        b.add(indices, _local_rhs);
    }

private:
    HCOpenBoundaryConditionData const& _data;
    std::size_t const _local_matrix_size;
};

}  // namespace ProcessLib
