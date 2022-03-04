/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Elements/MapBulkElementPoint.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

namespace ProcessLib
{
class SurfaceFluxLocalAssemblerInterface
{
public:
    virtual ~SurfaceFluxLocalAssemblerInterface() = default;

    virtual void integrate(
        std::size_t const element_id,
        std::vector<GlobalVector*> const& x,
        MeshLib::PropertyVector<double>& specific_flux,
        double const t,
        MeshLib::Mesh const& bulk_mesh,
        std::function<Eigen::Vector3d(std::size_t const,
                                      MathLib::Point3d const&, double const,
                                      std::vector<GlobalVector*> const&)>) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class SurfaceFluxLocalAssembler final
    : public SurfaceFluxLocalAssemblerInterface
{
protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

public:
    /// Precomputes the shape matrices for a given surface element.
    /// @param surface_element The surface element used for precomputing the
    /// shape matrices used later on for the integration.
    /// @param is_axially_symmetric corrects integration measure for cylinder
    /// coordinates.
    /// @param bulk_element_ids The id of the corresponding element in the bulk
    /// mesh.
    /// @param bulk_face_ids The id of the corresponding face in the bulk
    /// element.
    /// @param integration_order the order of the integration
    SurfaceFluxLocalAssembler(
        MeshLib::Element const& surface_element,
        std::size_t /*const local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        MeshLib::PropertyVector<std::size_t> const& bulk_element_ids,
        MeshLib::PropertyVector<std::size_t> const& bulk_face_ids)
        : _surface_element(surface_element),
          _integration_method(integration_order),
          _bulk_element_id(bulk_element_ids[surface_element.getID()]),
          _bulk_face_id(bulk_face_ids[surface_element.getID()])
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim, NumLib::ShapeMatrixType::N_J>(
                _surface_element, is_axially_symmetric, _integration_method);
        for (std::size_t ip = 0; ip < n_integration_points; ++ip)
        {
            _detJ_times_integralMeasure.push_back(
                shape_matrices[ip].detJ * shape_matrices[ip].integralMeasure);
        }
    }

    /// Integration for the element with the id \c element_id.
    /// @param element_id id of the element
    /// @param x The global vector containing the values for numerical
    /// integration.
    /// @param specific_flux PropertyVector the integration result will be
    /// stored into, where specific_flux has the same number of entries like
    /// mesh elements exists.
    /// @param t The integration is performed at the time t.
    /// @param bulk_mesh Stores a reference to the bulk mesh that is needed to
    /// fetch the information for the integration over the surface mesh.
    /// @param getFlux function that calculates the flux in the integration
    /// points of the face elements of the bulk element that belongs to the
    /// surface.
    void integrate(std::size_t const element_id,
                   std::vector<GlobalVector*> const& x,
                   MeshLib::PropertyVector<double>& specific_flux,
                   double const t,
                   MeshLib::Mesh const& bulk_mesh,
                   std::function<Eigen::Vector3d(
                       std::size_t const, MathLib::Point3d const&, double const,
                       std::vector<GlobalVector*> const&)>
                       getFlux) override
    {
        auto const& bulk_element = *bulk_mesh.getElement(_bulk_element_id);
        auto const surface_element_normal =
            getSurfaceNormal(_surface_element, bulk_element);

        double element_area = 0.0;
        std::size_t const n_integration_points =
            _integration_method.getNumberOfPoints();
        // specific_flux[id_of_element] +=
        //   int_{\Gamma_e} \alpha \cdot flux \cdot normal \d \Gamma /
        //   element_area
        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const& wp = _integration_method.getWeightedPoint(ip);

            auto const bulk_element_point =
                MeshLib::getBulkElementPoint(bulk_element, _bulk_face_id, wp);
            auto const bulk_flux =
                getFlux(_bulk_element_id, bulk_element_point, t, x);
            for (int component_id(0);
                 component_id < specific_flux.getNumberOfGlobalComponents();
                 ++component_id)
            {
                // TODO find solution for 2d case
                double const bulk_grad_times_normal(
                    Eigen::Map<Eigen::RowVectorXd const>(bulk_flux.data(),
                                                         bulk_flux.size())
                        .dot(surface_element_normal));

                specific_flux.getComponent(element_id, component_id) +=
                    bulk_grad_times_normal * _detJ_times_integralMeasure[ip] *
                    wp.getWeight();
            }
            element_area += _detJ_times_integralMeasure[ip] * wp.getWeight();
        }
        for (int component_id(0);
             component_id < specific_flux.getNumberOfGlobalComponents();
             ++component_id)
        {
            specific_flux.getComponent(element_id, component_id) /=
                element_area;
        }
    }

private:
    static Eigen::Vector3d getSurfaceNormal(
        MeshLib::Element const& surface_element,
        MeshLib::Element const& bulk_element)
    {
        Eigen::Vector3d surface_element_normal;

        if (surface_element.getGeomType() == MeshLib::MeshElemType::LINE)
        {
            auto const bulk_normal =
                MeshLib::FaceRule::getSurfaceNormal(bulk_element);
            auto const l0 = Eigen::Map<Eigen::Vector3d const>(
                surface_element.getNode(0)->getCoords());
            auto const l1 = Eigen::Map<Eigen::Vector3d const>(
                surface_element.getNode(1)->getCoords());
            Eigen::Vector3d const line = l1 - l0;
            surface_element_normal = line.cross(bulk_normal);
        }
        else
        {
            surface_element_normal =
                MeshLib::FaceRule::getSurfaceNormal(surface_element);
        }
        surface_element_normal.normalize();
        // At the moment (2016-09-28) the surface normal is not oriented
        // according to the right hand rule. Thus for an intuitive flux
        // output, i.e., inflow has positive sign, outflow has negative
        // sign, the normal must not be multiplied by -1.
        return surface_element_normal;
    };

    MeshLib::Element const& _surface_element;

    std::vector<double> _detJ_times_integralMeasure;

    IntegrationMethod const _integration_method;
    std::size_t const _bulk_element_id;
    std::size_t const _bulk_face_id;
};

}  // namespace ProcessLib
