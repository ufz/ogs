/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_CALCULATESURFACEFLUXLOCALASSEMBLER_H
#define PROCESSLIB_CALCULATESURFACEFLUXLOCALASSEMBLER_H

#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Parameter/Parameter.h"

#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "MapBulkElementPoint.h"

#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Elements/Elements.h"

namespace ProcessLib
{
class CalculateSurfaceFluxLocalAssemblerInterface
{
public:
    virtual ~CalculateSurfaceFluxLocalAssemblerInterface() = default;

    virtual void integrate(std::size_t element_id,
                           GlobalVector const& x,
                           MeshLib::PropertyVector<double>& balance,
                           Process const& bulk_process) = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class CalculateSurfaceFluxLocalAssembler final
    : public CalculateSurfaceFluxLocalAssemblerInterface
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
    CalculateSurfaceFluxLocalAssembler(
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
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &surface_element));

        std::size_t const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<typename ShapeMatricesType::ShapeMatrices> shape_matrices;
        shape_matrices.reserve(n_integration_points);
        _detJ_times_integralMeasure.reserve(n_integration_points);
        for (std::size_t ip = 0; ip < n_integration_points; ++ip)
        {
            shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                        ShapeFunction::NPOINTS);
            fe.template computeShapeFunctions<NumLib::ShapeMatrixType::N_J>(
                _integration_method.getWeightedPoint(ip).getCoords(),
                shape_matrices[ip], GlobalDim, is_axially_symmetric);
            _detJ_times_integralMeasure.push_back(
                shape_matrices[ip].detJ * shape_matrices[ip].integralMeasure);
        }
    }

    /// Integration for the element with the id \c element_id.
    /// @param element_id id of the element
    /// @param x The global vector containing the values for numerical
    /// integration.
    /// @param balance PropertyVector the integration result will be stored
    /// into, where balance has the same number of entries like mesh elements
    /// exists.
    /// @param bulk_process Reference to the bulk process is used to compute the
    /// flux.
    void integrate(std::size_t element_id,
                   GlobalVector const& x,
                   MeshLib::PropertyVector<double>& balance,
                   Process const& bulk_process) override
    {
        auto surface_element_normal =
            MeshLib::FaceRule::getSurfaceNormal(&_surface_element);
        surface_element_normal.normalize();
        // At the moment (2016-09-28) the surface normal is not oriented
        // according to the right hand rule
        // for correct results it is necessary to multiply the normal with -1
        surface_element_normal *= -1;

        std::size_t const n_integration_points =
            _integration_method.getNumberOfPoints();
        // balance[id_of_element] =
        //   int_{\Gamma_e} \alpha \cdot flux \cdot normal \d \Gamma
        for (std::size_t ip(0); ip < n_integration_points; ip++)
        {
            auto const& wp = _integration_method.getWeightedPoint(ip);

            auto const bulk_element_point = getBulkElementPoint(
                bulk_process.getMesh(), _bulk_element_id, _bulk_face_id, wp);
            auto const bulk_flux =
                bulk_process.getFlux(_bulk_element_id, bulk_element_point, x);

            for (std::size_t component_id(0);
                 component_id < balance.getNumberOfComponents();
                 ++component_id)
            {
                // TODO find solution for 2d case
                double const bulk_grad_times_normal(
                    Eigen::Map<Eigen::RowVectorXd const>(bulk_flux.data(),
                                                   bulk_flux.size())
                        .dot(Eigen::Map<Eigen::RowVectorXd const>(
                            surface_element_normal.getCoords(), 3)));

                balance.getComponent(element_id, component_id) +=
                    bulk_grad_times_normal * _detJ_times_integralMeasure[ip] *
                    wp.getWeight();
            }
        }
    }

private:
    MeshLib::Element const& _surface_element;

    std::vector<double> _detJ_times_integralMeasure;

    IntegrationMethod const _integration_method;
    std::size_t _bulk_element_id;
    std::size_t _bulk_face_id;
};

}  // ProcessLib

#endif  // PROCESSLIB_CALCULATESURFACEFLUXLOCALASSEMBLER_H
