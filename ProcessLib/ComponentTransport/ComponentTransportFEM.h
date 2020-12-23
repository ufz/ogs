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

#include <Eigen/Dense>
#include <numeric>
#include <vector>

#include "ComponentTransportProcessData.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/ProcessVariable.h"

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

    GlobalIndexType chemical_system_id = 0;

    double porosity = std::numeric_limits<double>::quiet_NaN();
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class ComponentTransportLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    ComponentTransportLocalAssemblerInterface() = default;

    void setStaggeredCoupledSolutions(
        std::size_t const /*mesh_item_id*/,
        CoupledSolutionsForStaggeredScheme* const coupling_term)
    {
        _coupled_solutions = coupling_term;
    }

    virtual void setChemicalSystemID(std::size_t const /*mesh_item_id*/) = 0;

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getInterpolatedLocalSolution(
        const double /*t*/,
        std::vector<GlobalVector*> const& int_pt_x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const = 0;

protected:
    CoupledSolutionsForStaggeredScheme* _coupled_solutions{nullptr};
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public ComponentTransportLocalAssemblerInterface
{
    // When monolithic scheme is adopted, nodal pressure and nodal concentration
    // are accessed by vector index.
    static const int pressure_index = 0;
    static const int first_concentration_index = ShapeFunction::NPOINTS;

    static const int pressure_size = ShapeFunction::NPOINTS;
    static const int concentration_size =
        ShapeFunction::NPOINTS;  // per component

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalBlockMatrixType =
        typename ShapeMatricesType::template MatrixType<pressure_size,
                                                        pressure_size>;
    using LocalSegmentVectorType =
        typename ShapeMatricesType::template VectorType<pressure_size>;

    using LocalMatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using LocalVectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    LocalAssemblerData(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ComponentTransportProcessData const& process_data,
        std::vector<std::reference_wrapper<ProcessVariable>>
            transport_process_variables)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _transport_process_variables(transport_process_variables)
    {
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(element, is_axially_symmetric,
                                                 _integration_method);
        auto const& medium =
            _process_data.media_map->getMedium(_element.getID());
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);

            pos.setIntegrationPoint(ip);

            _ip_data[ip].porosity =
                medium->property(MaterialPropertyLib::PropertyType::porosity)
                    .template initialValue<double>(pos, 0.);
        }
    }

    void setChemicalSystemID(std::size_t const /*mesh_item_id*/) override
    {
        assert(_process_data.chemical_solver_interface);
        // chemical system index map
        auto& chemical_system_index_map =
            _process_data.chemical_solver_interface->chemical_system_index_map;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].chemical_system_id = chemical_system_index_map.empty()
                                     ? 0
                                     : chemical_system_index_map.back() + 1;
            chemical_system_index_map.push_back(_ip_data[ip].chemical_system_id);
        }
    }

    void assemble(double const t, double const dt,
                  std::vector<double> const& local_x,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // Nodal DOFs include pressure
        int const num_nodal_dof = 1 + _transport_process_variables.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * num_nodal_dof);

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        // Get block matrices
        auto Kpp = local_K.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Mpp = local_M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Bp = local_b.template segment<pressure_size>(pressure_index);

        auto local_p = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);

        auto const number_of_components = num_nodal_dof - 1;
        for (int component_id = 0; component_id < number_of_components;
             ++component_id)
        {
            /*  Partitioned assembler matrix
             *  |  pp | pc1 | pc2 | pc3 |
             *  |-----|-----|-----|-----|
             *  | c1p | c1c1|  0  |  0  |
             *  |-----|-----|-----|-----|
             *  | c2p |  0  | c2c2|  0  |
             *  |-----|-----|-----|-----|
             *  | c3p |  0  |  0  | c3c3|
             */
            auto concentration_index =
                pressure_size + component_id * concentration_size;

            auto KCC =
                local_K.template block<concentration_size, concentration_size>(
                    concentration_index, concentration_index);
            auto MCC =
                local_M.template block<concentration_size, concentration_size>(
                    concentration_index, concentration_index);
            auto MCp =
                local_M.template block<concentration_size, pressure_size>(
                    concentration_index, pressure_index);
            auto MpC =
                local_M.template block<pressure_size, concentration_size>(
                    pressure_index, concentration_index);

            auto local_C = Eigen::Map<const NodalVectorType>(
                &local_x[concentration_index], concentration_size);

            assembleBlockMatrices(component_id, t, dt, local_C, local_p, KCC,
                                  MCC, MCp, MpC, Kpp, Mpp, Bp);
        }
    }

    void assembleBlockMatrices(
        int const component_id, double const t, double const dt,
        Eigen::Ref<const NodalVectorType> const& C_nodal_values,
        Eigen::Ref<const NodalVectorType> const& p_nodal_values,
        Eigen::Ref<LocalBlockMatrixType> KCC,
        Eigen::Ref<LocalBlockMatrixType> MCC,
        Eigen::Ref<LocalBlockMatrixType> MCp,
        Eigen::Ref<LocalBlockMatrixType> MpC,
        Eigen::Ref<LocalBlockMatrixType> Kpp,
        Eigen::Ref<LocalBlockMatrixType> Mpp,
        Eigen::Ref<LocalSegmentVectorType> Bp)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialPropertyLib::VariableArray vars;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        // Get material properties
        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        // Select the only valid for component transport liquid phase.
        auto const& phase = medium.phase("AqueousLiquid");

        // Assume that the component name is the same as the process variable
        // name. Components are shifted by one because the first one is always
        // pressure.
        auto const& component = phase.component(
            _transport_process_variables[component_id].get().getName());

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;
            auto& porosity = ip_data.porosity;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(C_nodal_values, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(p_nodal_values, N, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::concentration)] = C_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // update according to a particular porosity model
            porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t, dt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
                porosity;

            auto const& retardation_factor =
                component
                    .property(
                        MaterialPropertyLib::PropertyType::retardation_factor)
                    .template value<double>(vars, pos, t, dt);

            auto const& solute_dispersivity_transverse = medium.template value<
                double>(
                MaterialPropertyLib::PropertyType::transversal_dispersivity);

            auto const& solute_dispersivity_longitudinal =
                medium.template value<double>(
                    MaterialPropertyLib::PropertyType::
                        longitudinal_dispersivity);

            // Use the fluid density model to compute the density
            // TODO (renchao): concentration of which component as the argument
            // for calculation of fluid density
            auto const density =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);

            auto const decay_rate =
                component
                    .property(MaterialPropertyLib::PropertyType::decay_rate)
                    .template value<double>(vars, pos, t, dt);

            auto const& molecular_diffusion_coefficient =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    component
                        .property(MaterialPropertyLib::PropertyType::
                                      molecular_diffusion)
                        .value(vars, pos, t, dt));

            auto const& K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));

            // Use the viscosity model to compute the viscosity
            auto const mu =
                phase.property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);

            GlobalDimMatrixType const K_over_mu = K / mu;
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu *
                                          (dNdx * p_nodal_values - density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);

            const double drho_dp =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t, dt);

            const double drho_dC =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::concentration, pos,
                        t, dt);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const hydrodynamic_dispersion =
                velocity_magnitude != 0.0
                    ? GlobalDimMatrixType(porosity *
                                              molecular_diffusion_coefficient +
                                          solute_dispersivity_transverse *
                                              velocity_magnitude * I +
                                          (solute_dispersivity_longitudinal -
                                           solute_dispersivity_transverse) /
                                              velocity_magnitude * velocity *
                                              velocity.transpose())
                    : GlobalDimMatrixType(porosity *
                                              molecular_diffusion_coefficient +
                                          solute_dispersivity_transverse *
                                              velocity_magnitude * I);
            const double R_times_phi(retardation_factor * porosity);
            GlobalDimVectorType const mass_density_flow = velocity * density;
            auto const N_t_N = (N.transpose() * N).eval();
            if (_process_data.non_advective_form)
            {
                MCp.noalias() += N_t_N * (C_int_pt * R_times_phi * drho_dp * w);
                MCC.noalias() += N_t_N * (C_int_pt * R_times_phi * drho_dC * w);
                KCC.noalias() -= dNdx.transpose() * mass_density_flow * N * w;
            }
            else
            {
                KCC.noalias() +=
                    N.transpose() * mass_density_flow.transpose() * dNdx * w;
            }
            MCC.noalias() += N_t_N * (R_times_phi * density * w);
            KCC.noalias() += dNdx.transpose() * hydrodynamic_dispersion * dNdx *
                                 (density * w) +
                             N_t_N * (decay_rate * R_times_phi * density * w);

            MpC.noalias() += N_t_N * (porosity * drho_dC * w);

            // Calculate Mpp, Kpp, and bp in the first loop over components
            if (component_id == 0)
            {
                Mpp.noalias() += N_t_N * (porosity * drho_dp * w);
                Kpp.noalias() +=
                    dNdx.transpose() * K_over_mu * dNdx * (density * w);

                if (_process_data.has_gravity)
                {
                    Bp.noalias() += dNdx.transpose() * K_over_mu * b *
                                    (density * density * w);
                }
            }
        }
    }

    void assembleForStaggeredScheme(double const t, double const dt,
                                    Eigen::VectorXd const& local_x,
                                    Eigen::VectorXd const& local_xdot,
                                    int const process_id,
                                    std::vector<double>& local_M_data,
                                    std::vector<double>& local_K_data,
                                    std::vector<double>& local_b_data) override
    {
        if (process_id == _process_data.hydraulic_process_id)
        {
            assembleHydraulicEquation(t, dt, local_x, local_xdot, local_M_data,
                                      local_K_data, local_b_data);
        }
        else
        {
            // Go for assembling in an order of transport process id.
            assembleComponentTransportEquation(t, dt, local_x, local_xdot,
                                               local_M_data, local_K_data,
                                               local_b_data, process_id);
        }
    }

    void assembleHydraulicEquation(double const t,
                                   double const dt,
                                   Eigen::VectorXd const& local_x,
                                   Eigen::VectorXd const& local_xdot,
                                   std::vector<double>& local_M_data,
                                   std::vector<double>& local_K_data,
                                   std::vector<double>& local_b_data)
    {
        auto const local_p =
            local_x.template segment<pressure_size>(pressure_index);
        auto const local_C = local_x.template segment<concentration_size>(
            first_concentration_index);
        auto const local_Cdot =
            local_xdot.segment<concentration_size>(first_concentration_index);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, pressure_size, pressure_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, pressure_size, pressure_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, pressure_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");

        MaterialPropertyLib::VariableArray vars;

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;
            auto& porosity = ip_data.porosity;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_C, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::concentration)] = C_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // update according to a particular porosity model
            porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t, dt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
                porosity;

            // Use the fluid density model to compute the density
            // TODO: Concentration of which component as one of arguments for
            // calculation of fluid density
            auto const density =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);

            auto const& K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));

            // Use the viscosity model to compute the viscosity
            auto const mu =
                phase.property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);

            GlobalDimMatrixType const K_over_mu = K / mu;

            const double drho_dp =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::phase_pressure,
                        pos, t, dt);
            const double drho_dC =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template dValue<double>(
                        vars, MaterialPropertyLib::Variable::concentration, pos,
                        t, dt);

            // matrix assembly
            local_M.noalias() += w * N.transpose() * porosity * drho_dp * N;
            local_K.noalias() +=
                w * dNdx.transpose() * density * K_over_mu * dNdx;

            if (_process_data.has_gravity)
            {
                local_b.noalias() +=
                    w * density * density * dNdx.transpose() * K_over_mu * b;
            }

            // coupling term
            {
                double dot_C_int_pt = 0.0;
                NumLib::shapeFunctionInterpolate(local_Cdot, N, dot_C_int_pt);

                local_b.noalias() -=
                    w * N.transpose() * porosity * drho_dC * dot_C_int_pt;
            }
        }
    }

    void assembleComponentTransportEquation(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_xdot, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& /*local_b_data*/, int const transport_process_id)
    {
        auto const local_p =
            local_x.template segment<pressure_size>(pressure_index);
        auto const local_C = local_x.template segment<concentration_size>(
            first_concentration_index +
            (transport_process_id - 1) * concentration_size);
        auto const local_pdot =
            local_xdot.segment<pressure_size>(pressure_index);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, concentration_size, concentration_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, concentration_size, concentration_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialPropertyLib::VariableArray vars;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");
        // Hydraulic process id is 0 and thus transport process id starts
        // from 1.
        auto const component_id = transport_process_id - 1;
        auto const& component = phase.component(
            _transport_process_variables[component_id].get().getName());

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;
            auto& porosity = ip_data.porosity;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_C, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::concentration)] = C_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;

            // update according to a particular porosity model
            porosity =
                medium.property(MaterialPropertyLib::PropertyType::porosity)
                    .template value<double>(vars, pos, t, dt);
            vars[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
                porosity;

            auto const& retardation_factor =
                component
                    .property(
                        MaterialPropertyLib::PropertyType::retardation_factor)
                    .template value<double>(vars, pos, t, dt);

            auto const& solute_dispersivity_transverse = medium.template value<
                double>(
                MaterialPropertyLib::PropertyType::transversal_dispersivity);
            auto const& solute_dispersivity_longitudinal =
                medium.template value<double>(
                    MaterialPropertyLib::PropertyType::
                        longitudinal_dispersivity);

            // Use the fluid density model to compute the density
            auto const density =
                phase.property(MaterialPropertyLib::PropertyType::density)
                    .template value<double>(vars, pos, t, dt);
            auto const decay_rate =
                component
                    .property(MaterialPropertyLib::PropertyType::decay_rate)
                    .template value<double>(vars, pos, t, dt);

            auto const& molecular_diffusion_coefficient =
                MaterialPropertyLib::formEigenTensor<GlobalDim>(
                    component
                        .property(MaterialPropertyLib::PropertyType::
                                      molecular_diffusion)
                        .value(vars, pos, t, dt));

            auto const& K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));
            // Use the viscosity model to compute the viscosity
            auto const mu =
                phase.property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);

            GlobalDimMatrixType const K_over_mu = K / mu;
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu *
                                          (dNdx * local_p - density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * local_p);

            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const hydrodynamic_dispersion =
                velocity_magnitude != 0.0
                    ? GlobalDimMatrixType(porosity *
                                              molecular_diffusion_coefficient +
                                          solute_dispersivity_transverse *
                                              velocity_magnitude * I +
                                          (solute_dispersivity_longitudinal -
                                           solute_dispersivity_transverse) /
                                              velocity_magnitude * velocity *
                                              velocity.transpose())
                    : GlobalDimMatrixType(porosity *
                                              molecular_diffusion_coefficient +
                                          solute_dispersivity_transverse *
                                              velocity_magnitude * I);

            double const R_times_phi = retardation_factor * porosity;
            auto const N_t_N = (N.transpose() * N).eval();

            if (_process_data.non_advective_form)
            {
                const double drho_dC =
                    phase.property(MaterialPropertyLib::PropertyType::density)
                        .template dValue<double>(
                            vars, MaterialPropertyLib::Variable::concentration,
                            pos, t, dt);
                local_M.noalias() +=
                    N_t_N * (R_times_phi * C_int_pt * drho_dC * w);
            }

            local_M.noalias() += N_t_N * (R_times_phi * density * w);

            // coupling term

            if (_process_data.non_advective_form)
            {
                double dot_p_int_pt = 0.0;

                NumLib::shapeFunctionInterpolate(local_pdot, N, dot_p_int_pt);
                const double drho_dp =
                    phase.property(MaterialPropertyLib::PropertyType::density)
                        .template dValue<double>(
                            vars, MaterialPropertyLib::Variable::phase_pressure,
                            pos, t, dt);

                local_K.noalias() +=
                    N_t_N * ((R_times_phi * drho_dp * dot_p_int_pt) * w) -
                    dNdx.transpose() * velocity * N * (density * w);
            }
            else
            {
                local_K.noalias() +=
                    N.transpose() * velocity.transpose() * dNdx * (density * w);
            }
            local_K.noalias() +=
                dNdx.transpose() * hydrodynamic_dispersion * dNdx *
                    (density * w) +
                N_t_N * (decay_rate * R_times_phi * density * w);
        }
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const override
    {
        assert(x.size() == dof_table.size());

        auto const n_processes = x.size();
        std::vector<std::vector<double>> local_x;
        local_x.reserve(n_processes);

        for (std::size_t process_id = 0; process_id < n_processes; ++process_id)
        {
            auto const indices =
                NumLib::getIndices(_element.getID(), *dof_table[process_id]);
            assert(!indices.empty());
            local_x.push_back(x[process_id]->get(indices));
        }

        // only one process, must be monolithic.
        if (n_processes == 1)
        {
            auto const local_p = Eigen::Map<const NodalVectorType>(
                &local_x[0][pressure_index], pressure_size);
            auto const local_C = Eigen::Map<const NodalVectorType>(
                &local_x[0][first_concentration_index], concentration_size);
            return calculateIntPtDarcyVelocity(t, local_p, local_C, cache);
        }

        // multiple processes, must be staggered.
        {
            constexpr int pressure_process_id = 0;
            constexpr int concentration_process_id = 1;
            auto const local_p = Eigen::Map<const NodalVectorType>(
                &local_x[pressure_process_id][0], pressure_size);
            auto const local_C = Eigen::Map<const NodalVectorType>(
                &local_x[concentration_process_id][0], concentration_size);
            return calculateIntPtDarcyVelocity(t, local_p, local_C, cache);
        }
    }

    std::vector<double> const& calculateIntPtDarcyVelocity(
        const double t,
        Eigen::Ref<const NodalVectorType> const& p_nodal_values,
        Eigen::Ref<const NodalVectorType> const& C_nodal_values,
        std::vector<double>& cache) const
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialPropertyLib::VariableArray vars;

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& porosity = ip_data.porosity;

            pos.setIntegrationPoint(ip);

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(C_nodal_values, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(p_nodal_values, N, p_int_pt);

            vars[static_cast<int>(
                MaterialPropertyLib::Variable::concentration)] = C_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;
            vars[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
                porosity;

            // TODO (naumov) Temporary value not used by current material
            // models. Need extension of secondary variables interface.
            double const dt = std::numeric_limits<double>::quiet_NaN();
            auto const& K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
                medium.property(MaterialPropertyLib::PropertyType::permeability)
                    .value(vars, pos, t, dt));
            auto const mu =
                phase.property(MaterialPropertyLib::PropertyType::viscosity)
                    .template value<double>(vars, pos, t, dt);
            GlobalDimMatrixType const K_over_mu = K / mu;

            cache_mat.col(ip).noalias() = -K_over_mu * dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                auto const rho_w =
                    phase.property(MaterialPropertyLib::PropertyType::density)
                        .template value<double>(vars, pos, t, dt);
                auto const b = _process_data.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                cache_mat.col(ip).noalias() += K_over_mu * rho_w * b;
            }
        }

        return cache;
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override
    {
        auto const local_p = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);
        auto const local_C = Eigen::Map<const NodalVectorType>(
            &local_x[first_concentration_index], concentration_size);

        // Eval shape matrices at given point
        // Note: Axial symmetry is set to false here, because we only need dNdx
        // here, which is not affected by axial symmetry.
        auto const shape_matrices =
            NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                         GlobalDim>(
                _element, false /*is_axially_symmetric*/,
                std::array{pnt_local_coords})[0];

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialPropertyLib::VariableArray vars;

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");

        // local_x contains the local concentration and pressure values
        double c_int_pt;
        NumLib::shapeFunctionInterpolate(local_C, shape_matrices.N, c_int_pt);
        vars[static_cast<int>(MaterialPropertyLib::Variable::concentration)]
            .emplace<double>(c_int_pt);

        double p_int_pt;
        NumLib::shapeFunctionInterpolate(local_p, shape_matrices.N, p_int_pt);
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)]
            .emplace<double>(p_int_pt);

        // TODO (naumov) Temporary value not used by current material models.
        // Need extension of secondary variables interface.
        double const dt = std::numeric_limits<double>::quiet_NaN();
        auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium.property(MaterialPropertyLib::PropertyType::permeability)
                .value(vars, pos, t, dt));

        auto const mu =
            phase.property(MaterialPropertyLib::PropertyType::viscosity)
                .template value<double>(vars, pos, t, dt);
        GlobalDimMatrixType const K_over_mu = K / mu;

        GlobalDimVectorType q = -K_over_mu * shape_matrices.dNdx * local_p;
        auto const rho_w =
            phase.property(MaterialPropertyLib::PropertyType::density)
                .template value<double>(vars, pos, t, dt);
        if (_process_data.has_gravity)
        {
            auto const b = _process_data.specific_body_force;
            q += K_over_mu * rho_w * b;
        }
        Eigen::Vector3d flux(0.0, 0.0, 0.0);
        flux.head<GlobalDim>() = rho_w * q;
        return flux;
    }

    std::vector<double> interpolateNodalValuesToIntegrationPoints(
        std::vector<double> const& local_x) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> interpolated_values(n_integration_points);
        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            NumLib::shapeFunctionInterpolate(local_x, _ip_data[ip].N,
                                             interpolated_values[ip]);
        }
        return interpolated_values;
    }

    std::vector<double> const& getInterpolatedLocalSolution(
        const double /*t*/,
        std::vector<GlobalVector*> const& int_pt_x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        assert(_process_data.chemical_process_data);
        assert(int_pt_x.size() == 1);

        cache.clear();
        auto const ele_id = _element.getID();
        auto const& indices = _process_data.chemical_process_data
                                  ->chemical_system_index_map[ele_id];
        cache = int_pt_x[0]->get(indices);

        return cache;
    }

    void computeSecondaryVariableConcrete(
        double const t,
        double const /*dt*/,
        Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& /*local_x_dot*/) override
    {
        auto const local_p =
            local_x.template segment<pressure_size>(pressure_index);
        auto const local_C = local_x.template segment<concentration_size>(
            first_concentration_index);

        std::vector<double> ele_velocity;
        calculateIntPtDarcyVelocity(t, local_p, local_C, ele_velocity);

        auto const n_integration_points =
            _integration_method.getNumberOfPoints();
        auto const ele_velocity_mat =
            MathLib::toMatrix(ele_velocity, GlobalDim, n_integration_points);

        auto ele_id = _element.getID();
        Eigen::Map<LocalVectorType>(
            &(*_process_data.mesh_prop_velocity)[ele_id * GlobalDim],
            GlobalDim) =
            ele_velocity_mat.rowwise().sum() / n_integration_points;
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              double const /*t*/, double const /*dt*/) override
    {
    }

private:
    MeshLib::Element const& _element;
    ComponentTransportProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<std::reference_wrapper<ProcessVariable>> const
        _transport_process_variables;

    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
