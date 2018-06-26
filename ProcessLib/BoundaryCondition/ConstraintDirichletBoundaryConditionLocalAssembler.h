/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "NumLib/DOF/DOFTableUtil.h"

#include "ProcessLib/Process.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "MeshLib/Elements/MapBulkElementPoint.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Elements/Elements.h"

namespace ProcessLib
{
struct IntegrationPointData final
{
    IntegrationPointData(double const& detJ,
                         double const& integral_measure,
                         double const& integration_weight)
        : detJ_times_integralMeasure_times_weight(detJ * integral_measure *
                                                  integration_weight)
    {
    }

    double const detJ_times_integralMeasure_times_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class ConstraintDirichletBoundaryConditionLocalAssemblerInterface
{
public:
    virtual ~ConstraintDirichletBoundaryConditionLocalAssemblerInterface() =
        default;

    virtual double integrate(
        GlobalVector const& x, double const t, MeshLib::Mesh const& bulk_mesh,
        std::function<Eigen::Vector3d(std::size_t const,
                                      MathLib::Point3d const&, double const,
                                      GlobalVector const&)> const& getFlux) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class ConstraintDirichletBoundaryConditionLocalAssembler final
    : public ConstraintDirichletBoundaryConditionLocalAssemblerInterface
{
protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

public:
    /// Precomputes the shape matrices for a given surface element.
    /// @param surface_element The surface element used for precomputing the
    /// @param local_matrix_size Unused parameter, only necessary for the
    /// creation scheme of local assemblers
    /// shape matrices used later on for the integration.
    /// @param is_axially_symmetric Corrects integration measure for cylinder
    /// coordinates.
    /// @param integration_order The order of the integration.
    /// @param bulk_ids Pairs of bulk element ids and bulk element face ids.
    ConstraintDirichletBoundaryConditionLocalAssembler(
        MeshLib::Element const& surface_element,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric, unsigned const integration_order,
        std::vector<std::pair<std::size_t, unsigned>> bulk_ids)
        : _surface_element(surface_element),
          _integration_method(integration_order),
          _bulk_element_id(bulk_ids[_surface_element.getID()].first),
          _bulk_face_id(bulk_ids[_surface_element.getID()].second)
    {
        (void)local_matrix_size; // unused, but needed for the interface
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &_surface_element));

        std::size_t const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<
            typename ShapeMatricesType::ShapeMatrices,
            Eigen::aligned_allocator<typename ShapeMatricesType::ShapeMatrices>>
            shape_matrices;
        shape_matrices.reserve(n_integration_points);
        _ip_data.reserve(n_integration_points);
        for (std::size_t ip = 0; ip < n_integration_points; ++ip)
        {
            shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                        ShapeFunction::NPOINTS);
            fe.template computeShapeFunctions<NumLib::ShapeMatrixType::N_J>(
                _integration_method.getWeightedPoint(ip).getCoords(),
                shape_matrices[ip], GlobalDim, is_axially_symmetric);
            auto const& wp = _integration_method.getWeightedPoint(ip);
            _ip_data.emplace_back(shape_matrices[ip].detJ,
                                  shape_matrices[ip].integralMeasure,
                                  wp.getWeight());
        }
    }

    /// Integration for the element with the id \c element_id.
    /// @param x The global vector containing the values for numerical
    /// integration.
    /// @param t The point in time the the integration will be performed.
    /// @param bulk_mesh The bulk mesh the process is defined on.
    /// @param getFlux The function of the constraining process used to
    /// calculate the flux.
    double integrate(
        GlobalVector const& x, double const t, MeshLib::Mesh const& bulk_mesh,
        std::function<Eigen::Vector3d(
            std::size_t const, MathLib::Point3d const&, double const,
            GlobalVector const&)> const& getFlux) override
    {
        MathLib::Vector3 surface_element_normal;
        if (_surface_element.getDimension() < 2)
        {
            auto const bulk_element_normal =
                MeshLib::FaceRule::getSurfaceNormal(
                    bulk_mesh.getElement(_bulk_element_id));
            MathLib::Vector3 const edge_vector(*_surface_element.getNode(0),
                                               *_surface_element.getNode(1));
            surface_element_normal =
                MathLib::crossProduct(bulk_element_normal, edge_vector);
        } else {
            surface_element_normal =
                MeshLib::FaceRule::getSurfaceNormal(&_surface_element);
        }

        surface_element_normal.normalize();
        // At the moment (2018-04-26) the surface normal is not oriented
        // according to the right hand rule
        // for correct results it is necessary to multiply the normal with -1
        surface_element_normal *= -1;

        auto const n_integration_points =
            _integration_method.getNumberOfPoints();
        // integrated_value +=
        //   int_{\Gamma_e} \alpha \cdot flux \cdot normal \d \Gamma
        double integrated_value = 0;
        for (auto ip(0); ip < n_integration_points; ip++)
        {
            auto const& wp = _integration_method.getWeightedPoint(ip);
            auto const bulk_element_point = MeshLib::getBulkElementPoint(
                bulk_mesh, _bulk_element_id, _bulk_face_id, wp);
            auto const bulk_flux =
                getFlux(_bulk_element_id, bulk_element_point, t, x);

            // TODO find solution for 2d case
            double const bulk_grad_times_normal(
                Eigen::Map<Eigen::RowVectorXd const>(bulk_flux.data(),
                                                     bulk_flux.size())
                    .dot(Eigen::Map<Eigen::RowVectorXd const>(
                        surface_element_normal.getCoords(), 3)));

            integrated_value +=
                bulk_grad_times_normal *
                _ip_data[ip].detJ_times_integralMeasure_times_weight;
        }
        return integrated_value;
    }

private:
    MeshLib::Element const& _surface_element;

    std::vector<IntegrationPointData> _ip_data;

    IntegrationMethod const _integration_method;
    std::size_t const _bulk_element_id;
    unsigned const _bulk_face_id;
};

}  // ProcessLib
