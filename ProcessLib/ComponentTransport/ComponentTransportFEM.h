/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>


#include "ComponentTransportProcessData.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace ComponentTransport
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         GlobalDimNodalMatrixType const& dNdx_,
                         double const& integration_weight_)
        : N(N_), dNdx(dNdx_), integration_weight(integration_weight_)
    {}

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;
};

const unsigned NUM_NODAL_DOF = 2;

class ComponentTransportLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    virtual std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const = 0;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public ComponentTransportLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalMatrixType = typename ShapeMatricesType::template MatrixType<
        NUM_NODAL_DOF * ShapeFunction::NPOINTS,
        NUM_NODAL_DOF * ShapeFunction::NPOINTS>;
    using LocalVectorType =
        typename ShapeMatricesType::template VectorType<NUM_NODAL_DOF *
                                                        ShapeFunction::NPOINTS>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    LocalAssemblerData(MeshLib::Element const& element,
                       std::size_t const local_matrix_size,
                       bool is_axially_symmetric,
                       unsigned const integration_order,
                       ComponentTransportProcessData const& process_data)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _darcy_velocities(
              GlobalDim,
              std::vector<double>(_integration_method.getNumberOfPoints()))
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const num_nodes = ShapeFunction::NPOINTS;
        auto p_nodal_values =
            Eigen::Map<const NodalVectorType>(&local_x[num_nodes], num_nodes);

        auto const& b = _process_data.specific_body_force;

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        auto KCC = local_K.template block<num_nodes, num_nodes>(0, 0);
        auto MCC = local_M.template block<num_nodes, num_nodes>(0, 0);
        auto Kpp =
            local_K.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Mpp =
            local_M.template block<num_nodes, num_nodes>(num_nodes, num_nodes);
        auto Bp = local_b.template block<num_nodes, 1>(num_nodes, 0);

        for (std::size_t ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            // \todo the argument to getValue() has to be changed for non
            // constant storage model
            auto const specific_storage =
                _process_data.porous_media_properties.getSpecificStorage(t, pos)
                    .getValue(0.0);

            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;
            // Order matters: First C, then p!
            NumLib::shapeFunctionInterpolate(local_x, N, C_int_pt, p_int_pt);

            // \todo the first argument has to be changed for non constant
            // porosity model
            auto const porosity =
                _process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(0.0, C_int_pt);

            auto const retardation_factor =
                _process_data.retardation_factor(t, pos)[0];

            auto const& solute_dispersivity_transverse =
                _process_data.solute_dispersivity_transverse(t, pos)[0];
            auto const& solute_dispersivity_longitudinal =
                _process_data.solute_dispersivity_longitudinal(t, pos)[0];

            // Use the fluid density model to compute the density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density_water = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);
            auto const& decay_rate = _process_data.decay_rate(t, pos)[0];
            auto const& molecular_diffusion_coefficient =
                _process_data.molecular_diffusion_coefficient(t, pos)[0];

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos);
            // Use the viscosity model to compute the viscosity
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

            GlobalDimMatrixType const K_over_mu = K / mu;

            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                    ? GlobalDimVectorType(
                          -perm_over_visc *
                          (dNdx * p_nodal_values - density_water * b))
                    : GlobalDimVectorType(-perm_over_visc * dNdx *
                                          p_nodal_values);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const hydrodynamic_dispersion =
                velocity_magnitude != 0.0
                    ? GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                              I +
                          (solute_dispersivity_longitudinal -
                           solute_dispersivity_transverse) /
                              velocity_magnitude * velocity *
                              velocity.transpose())
                    : GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                          I);

            // matrix assembly
            KCC.noalias() +=
                (dNdx.transpose() * hydrodynamic_dispersion * dNdx +
                 N.transpose() * velocity.transpose() * dNdx +
                 N.transpose() * decay_rate * porosity * retardation_factor *
                     N) *
                w;
            MCC.noalias() +=
                w * N.transpose() * porosity * retardation_factor * N;
            Kpp.noalias() += w * dNdx.transpose() * perm_over_visc * dNdx;
            Mpp.noalias() += w * N.transpose() * specific_storage * N;
            if (_process_data.has_gravity)
                Bp += w * density_water * dNdx.transpose() * perm_over_visc * b;
            /* with Oberbeck-Boussing assumption density difference only exists
             * in buoyancy effects */
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, std::vector<double> const& local_x) override
    {
        SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& K =
            _process_data.porous_media_properties.getIntrinsicPermeability(t,
                                                                           pos);
        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        auto const mu = _process_data.fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
        GlobalDimMatrixType const K_over_mu = K / mu;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[ShapeFunction::NPOINTS], ShapeFunction::NPOINTS);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;

            GlobalDimVectorType velocity = -K_over_mu * dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                double C_int_pt = 0.0;
                double p_int_pt = 0.0;
                NumLib::shapeFunctionInterpolate(local_x, N, C_int_pt,
                                                 p_int_pt);
                vars[static_cast<int>(
                    MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
                vars[static_cast<int>(
                    MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;

                auto const rho_w = _process_data.fluid_properties->getValue(
                    MaterialLib::Fluid::FluidPropertyType::Density, vars);
                auto const b = _process_data.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                velocity += K_over_mu * rho_w * b;
            }

            for (unsigned d = 0; d < GlobalDim; ++d)
            {
                _darcy_velocities[d][ip] = velocity[d];
            }
        }
    }


    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getIntPtDarcyVelocityX(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 0);
        return _darcy_velocities[0];
    }

    std::vector<double> const& getIntPtDarcyVelocityY(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 1);
        return _darcy_velocities[1];
    }

    std::vector<double> const& getIntPtDarcyVelocityZ(
        std::vector<double>& /*cache*/) const override
    {
        assert(_darcy_velocities.size() > 2);
        return _darcy_velocities[2];
    }

private:
    MeshLib::Element const& _element;
    ComponentTransportProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;
    std::vector<std::vector<double>> _darcy_velocities;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
