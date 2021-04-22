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

#include <Eigen/Dense>
#include <vector>

#include "HTLocalAssemblerInterface.h"
#include "HTProcessData.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib
{
namespace HT
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class HTFEM : public HTLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    HTFEM(MeshLib::Element const& element,
          std::size_t const local_matrix_size,
          bool const is_axially_symmetric,
          unsigned const integration_order,
          HTProcessData const& process_data,
          const unsigned dof_per_node)
        : HTLocalAssemblerInterface(),
          _element(element),
          _process_data(process_data),
          _integration_method(integration_order)
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * dof_per_node);
        (void)local_matrix_size;
        (void)dof_per_node;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    /// Computes the flux in the point \c pnt_local_coords that is given in
    /// local coordinates using the values from \c local_x.
    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override
    {
        // Eval shape matrices at given point
        // Note: Axial symmetry is set to false here, because we only need dNdx
        // here, which is not affected by axial symmetry.
        auto const shape_matrices =
            NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                         GlobalDim>(
                _element, false /*is_axially_symmetric*/,
                std::array{pnt_local_coords})[0];

        ParameterLib::SpatialPosition pos;
        pos.setElementID(this->_element.getID());

        MaterialPropertyLib::VariableArray vars;

        // local_x contains the local temperature and pressure values
        double T_int_pt = 0.0;
        double p_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, shape_matrices.N, T_int_pt,
                                         p_int_pt);

        vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
            T_int_pt;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            p_int_pt;

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");

        // TODO (naumov) Temporary value not used by current material models.
        // Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();
        // fetch permeability, viscosity, density
        auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium.property(MaterialPropertyLib::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const mu =
            liquid_phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);
        GlobalDimMatrixType const K_over_mu = K / mu;

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_x[local_x.size() / 2], ShapeFunction::NPOINTS);
        GlobalDimVectorType q =
            -K_over_mu * shape_matrices.dNdx * p_nodal_values;

        if (this->_process_data.has_gravity)
        {
            auto const rho_w =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);
            auto const b = this->_process_data.specific_body_force;
            q += K_over_mu * rho_w * b;
        }

        Eigen::Vector3d flux;
        flux.head<GlobalDim>() = q;
        return flux;
    }

protected:
    MeshLib::Element const& _element;
    HTProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;

    double getHeatEnergyCoefficient(
        MaterialPropertyLib::VariableArray const& vars, const double porosity,
        const double fluid_density, const double specific_heat_capacity_fluid,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt)
    {
        auto const& medium =
            *_process_data.media_map->getMedium(this->_element.getID());
        auto const& solid_phase = medium.phase("Solid");

        auto const specific_heat_capacity_solid =
            solid_phase
                .property(
                    MaterialPropertyLib::PropertyType::specific_heat_capacity)
                .template value<double>(vars, pos, t, dt);

        auto const solid_density =
            solid_phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);

        return solid_density * specific_heat_capacity_solid * (1 - porosity) +
               fluid_density * specific_heat_capacity_fluid * porosity;
    }

    GlobalDimMatrixType getThermalConductivityDispersivity(
        MaterialPropertyLib::VariableArray const& vars,
        const double fluid_density, const double specific_heat_capacity_fluid,
        const GlobalDimVectorType& velocity, const GlobalDimMatrixType& I,
        ParameterLib::SpatialPosition const& pos, double const t,
        double const dt)
    {
        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());

        auto const thermal_conductivity =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium
                    .property(
                        MaterialPropertyLib::PropertyType::thermal_conductivity)
                    .value(vars, pos, t, dt));

        auto const thermal_dispersivity_longitudinal =
            medium
                .property(MaterialPropertyLib::PropertyType::
                              thermal_longitudinal_dispersivity)
                .template value<double>();
        auto const thermal_dispersivity_transversal =
            medium
                .property(MaterialPropertyLib::PropertyType::
                              thermal_transversal_dispersivity)
                .template value<double>();

        double const velocity_magnitude = velocity.norm();

        if (velocity_magnitude < std::numeric_limits<double>::epsilon())
        {
            return thermal_conductivity;
        }

        GlobalDimMatrixType const thermal_dispersivity =
            fluid_density * specific_heat_capacity_fluid *
            (thermal_dispersivity_transversal * velocity_magnitude * I +
             (thermal_dispersivity_longitudinal -
              thermal_dispersivity_transversal) /
                 velocity_magnitude * velocity * velocity.transpose());

        return thermal_conductivity + thermal_dispersivity;
    }

    std::vector<double> const& getIntPtDarcyVelocityLocal(
        const double t, std::vector<double> const& local_x,
        std::vector<double>& cache) const
    {
        std::vector<double> local_p{
            local_x.data() + pressure_index,
            local_x.data() + pressure_index + pressure_size};
        std::vector<double> local_T{
            local_x.data() + temperature_index,
            local_x.data() + temperature_index + temperature_size};

        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialPropertyLib::VariableArray vars;

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_p[0], ShapeFunction::NPOINTS);

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& liquid_phase = medium.phase("AqueousLiquid");

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;

            pos.setIntegrationPoint(ip);

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);
            NumLib::shapeFunctionInterpolate(local_T, N, T_int_pt);

            vars[static_cast<int>(MaterialPropertyLib::Variable::temperature)] =
                T_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // TODO (naumov) Temporary value not used by current material
            // models. Need extension of secondary variables interface.
            double const dt = std::numeric_limits<double>::quiet_NaN();
            auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));

            auto const mu =
                liquid_phase
                    .property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);
            GlobalDimMatrixType const K_over_mu = K / mu;

            cache_mat.col(ip).noalias() = -K_over_mu * dNdx * p_nodal_values;

            if (_process_data.has_gravity)
            {
                auto const rho_w =
                    liquid_phase
                        .property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, t, dt);
                auto const b = _process_data.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                cache_mat.col(ip).noalias() += K_over_mu * rho_w * b;
            }
        }

        return cache;
    }

protected:
    static const int pressure_index = ShapeFunction::NPOINTS;
    static const int pressure_size = ShapeFunction::NPOINTS;
    static const int temperature_index = 0;
    static const int temperature_size = ShapeFunction::NPOINTS;
};

}  // namespace HT
}  // namespace ProcessLib
