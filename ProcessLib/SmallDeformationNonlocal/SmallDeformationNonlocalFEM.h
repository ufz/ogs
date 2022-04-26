/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>
#include <limits>
#include <memory>
#include <vector>

#include "Damage.h"
#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MaterialLib/SolidModels/Ehlers.h"
#include "MaterialLib/SolidModels/SelectSolidConstitutiveRelation.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "MeshLib/findElementsWithinRadius.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/Deformation/Divergence.h"
#include "ProcessLib/Deformation/LinearBMatrix.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/SetOrGetIntegrationPointData.h"
#include "SmallDeformationNonlocalProcessData.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
namespace MPL = MaterialPropertyLib;

/// Used for the extrapolation of the integration point values. It is ordered
/// (and stored) by integration points.
template <typename ShapeMatrixType>
struct SecondaryData
{
    std::vector<ShapeMatrixType, Eigen::aligned_allocator<ShapeMatrixType>> N;
};

template <typename ShapeFunction, typename IntegrationMethod,
          int DisplacementDim>
class SmallDeformationNonlocalLocalAssembler
    : public SmallDeformationNonlocalLocalAssemblerInterface<DisplacementDim>
{
public:
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;
    using BMatricesType = BMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using BMatrixType = typename BMatricesType::BMatrixType;
    using StiffnessMatrixType = typename BMatricesType::StiffnessMatrixType;
    using NodalForceVectorType = typename BMatricesType::NodalForceVectorType;
    using NodalDisplacementVectorType =
        typename BMatricesType::NodalForceVectorType;
    using IpData =
        IntegrationPointData<BMatricesType, ShapeMatricesType, DisplacementDim>;

    SmallDeformationNonlocalLocalAssembler(
        SmallDeformationNonlocalLocalAssembler const&) = delete;
    SmallDeformationNonlocalLocalAssembler(
        SmallDeformationNonlocalLocalAssembler&&) = delete;

    SmallDeformationNonlocalLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        SmallDeformationNonlocalProcessData<DisplacementDim>& process_data)
        : _process_data(process_data),
          _integration_method(integration_order),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        _ip_data.reserve(n_integration_points);
        _secondary_data.N.resize(n_integration_points);

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      DisplacementDim>(e, is_axially_symmetric,
                                                       _integration_method);

        auto& solid_material =
            MaterialLib::Solids::selectSolidConstitutiveRelation(
                _process_data.solid_materials,
                _process_data.material_ids,
                e.getID());
        auto* ehlers_solid_material = dynamic_cast<
            MaterialLib::Solids::Ehlers::SolidEhlers<DisplacementDim>*>(
            &solid_material);
        if (ehlers_solid_material == nullptr)
        {
            OGS_FATAL(
                "The SmallDeformationNonlocalLocal process supports only "
                "Ehlers material at the moment. For other materials the "
                "interface must be extended first.");
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(*ehlers_solid_material);
            auto& ip_data = _ip_data[ip];
            auto const& sm = shape_matrices[ip];
            _ip_data[ip].integration_weight =
                _integration_method.getWeightedPoint(ip).getWeight() *
                sm.integralMeasure * sm.detJ;

            ip_data.N = sm.N;
            ip_data.dNdx = sm.dNdx;

            // Initialize current time step values
            ip_data.sigma.setZero(
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim));
            ip_data.eps.setZero(MathLib::KelvinVector::kelvin_vector_dimensions(
                DisplacementDim));

            // Previous time step values are not initialized and are set later.
            ip_data.sigma_prev.resize(
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim));
            ip_data.eps_prev.resize(
                MathLib::KelvinVector::kelvin_vector_dimensions(
                    DisplacementDim));

            _secondary_data.N[ip] = shape_matrices[ip].N;

            ip_data.coordinates = getSingleIntegrationPointCoordinates(ip);
        }
    }

    std::size_t setIPDataInitialConditions(std::string const& name,
                                    double const* values,
                                    int const integration_order) override
    {
        if (integration_order !=
            static_cast<int>(_integration_method.getIntegrationOrder()))
        {
            OGS_FATAL(
                "Setting integration point initial conditions; The integration "
                "order of the local assembler for element {:d} is different "
                "from the integration order in the initial condition.",
                _element.getID());
        }

        if (name == "sigma_ip")
        {
            return setSigma(values);
        }

        if (name == "kappa_d_ip")
        {
            return ProcessLib::setIntegrationPointScalarData(values, _ip_data,
                                                             &IpData::kappa_d);
        }

        return 0;
    }

    void setIPDataInitialConditionsFromCellData(
        std::string const& name, std::vector<double> const& value) override
    {
        if (name == "kappa_d_ip")
        {
            if (value.size() != 1)
            {
                OGS_FATAL(
                    "CellData for kappa_d initial conditions has wrong number "
                    "of components. 1 expected, got {:d}.",
                    value.size());
            }
            setKappaD(value[0]);
        }
    }

    double alpha_0(double const distance2) const
    {
        double const internal_length2 = _process_data.internal_length_squared;
        return (distance2 > internal_length2)
                   ? 0
                   : (1 - distance2 / (internal_length2)) *
                         (1 - distance2 / (internal_length2));
    }

    void nonlocal(
        std::size_t const /*mesh_item_id*/,
        std::vector<
            std::unique_ptr<SmallDeformationNonlocalLocalAssemblerInterface<
                DisplacementDim>>> const& local_assemblers) override
    {
        auto const search_element_ids = MeshLib::findElementsWithinRadius(
            _element, _process_data.internal_length_squared);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> distances;  // Cache for ip-ip distances.
        //
        // For every integration point in this element collect the neighbouring
        // integration points falling in given radius (internal length) and
        // compute the alpha_kl weight.
        //
        for (unsigned k = 0; k < n_integration_points; k++)
        {
            //
            // Collect the integration points.
            //

            auto const& xyz = _ip_data[k].coordinates;

            // For all neighbors of element
            for (auto const search_element_id : search_element_ids)
            {
                auto const& la = local_assemblers[search_element_id];
                la->getIntegrationPointCoordinates(xyz, distances);
                for (int ip = 0; ip < static_cast<int>(distances.size()); ++ip)
                {
                    if (distances[ip] >= _process_data.internal_length_squared)
                    {
                        continue;
                    }
                    // save into current ip_k
                    _ip_data[k].non_local_assemblers.push_back(
                        {la->getIPDataPtr(ip),
                         std::numeric_limits<double>::quiet_NaN(),
                         distances[ip]});
                }
            }
            if (_ip_data[k].non_local_assemblers.size() == 0)
            {
                OGS_FATAL("no neighbours found!");
            }

            double a_k_sum_m = 0;
            for (auto const& tuple : _ip_data[k].non_local_assemblers)
            {
                double const distance2_m = tuple.distance2;

                auto const& w_m = tuple.ip_l_pointer->integration_weight;

                a_k_sum_m += w_m * alpha_0(distance2_m);
            }

            //
            // Calculate alpha_kl =
            //       alpha_0(|x_k - x_l|) / int_{m \in ip} alpha_0(|x_k - x_m|)
            //
            for (auto& tuple : _ip_data[k].non_local_assemblers)
            {
                double const distance2_l = tuple.distance2;
                double const a_kl = alpha_0(distance2_l) / a_k_sum_m;

                // Store the a_kl already multiplied with the integration
                // weight of that l integration point.
                auto const w_l = tuple.ip_l_pointer->integration_weight;
                tuple.alpha_kl_times_w_l = a_kl * w_l;
            }
        }
    }

    Eigen::Vector3d getSingleIntegrationPointCoordinates(
        int integration_point) const
    {
        auto const& N = _secondary_data.N[integration_point];

        Eigen::Vector3d xyz = Eigen::Vector3d::Zero();  // Resulting coordinates
        auto* nodes = _element.getNodes();
        for (int i = 0; i < N.size(); ++i)
        {
            auto const& node_coordinates{nodes[i]->asEigenVector3d()};
            xyz += node_coordinates * N[i];
        }
        return xyz;
    }

    /// For each of the current element's integration points the squared
    /// distance from the current integration point is computed and stored in
    /// the given distances cache.
    void getIntegrationPointCoordinates(
        Eigen::Vector3d const& coords,
        std::vector<double>& distances) const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        distances.resize(n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& xyz = _ip_data[ip].coordinates;
            distances[ip] = (xyz - coords).squaredNorm();
        }
    }

    void assemble(double const /*t*/, double const /*dt*/,
                  std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_xdot*/,
                  std::vector<double>& /*local_M_data*/,
                  std::vector<double>& /*local_K_data*/,
                  std::vector<double>& /*local_b_data*/) override
    {
        OGS_FATAL(
            "SmallDeformationNonlocalLocalAssembler: assembly without jacobian "
            "is not "
            "implemented.");
    }

    void preAssemble(double const t, double const dt,
                     std::vector<double> const& local_x) override
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        MPL::VariableArray variables;
        MPL::VariableArray variables_prev;
        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(_element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);
            auto const& eps_prev = _ip_data[ip].eps_prev;
            auto const& sigma_prev = _ip_data[ip].sigma_prev;

            auto& eps = _ip_data[ip].eps;
            auto& sigma = _ip_data[ip].sigma;
            auto& C = _ip_data[ip].C;
            auto& state = _ip_data[ip].material_state_variables;
            double const& damage_prev = _ip_data[ip].damage_prev;

            eps.noalias() =
                B *
                Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
                    local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

            // sigma is for plastic part only.
            std::unique_ptr<
                MathLib::KelvinVector::KelvinMatrixType<DisplacementDim>>
                new_C;
            std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
                DisplacementDim>::MaterialStateVariables>
                new_state;

            // Compute sigma_eff from damage total stress sigma
            using KelvinVectorType = typename BMatricesType::KelvinVectorType;
            KelvinVectorType const sigma_eff_prev =
                sigma_prev /
                (1. - damage_prev);  // damage_prev is in [0,1) range. See
                                     // calculateDamage() function.

            variables_prev[static_cast<int>(MPL::Variable::stress)]
                .emplace<
                    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                    sigma_eff_prev);
            variables_prev[static_cast<int>(MPL::Variable::mechanical_strain)]
                .emplace<
                    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                    eps_prev);
            variables_prev[static_cast<int>(MPL::Variable::temperature)]
                .emplace<double>(_process_data.reference_temperature);
            variables[static_cast<int>(MPL::Variable::mechanical_strain)]
                .emplace<
                    MathLib::KelvinVector::KelvinVectorType<DisplacementDim>>(
                    eps);
            variables[static_cast<int>(MPL::Variable::temperature)]
                .emplace<double>(_process_data.reference_temperature);

            auto&& solution = _ip_data[ip].solid_material.integrateStress(
                variables_prev, variables, t, x_position, dt, *state);

            if (!solution)
            {
                OGS_FATAL("Computation of local constitutive relation failed.");
            }

            std::tie(sigma, state, C) = std::move(*solution);

            /// Compute only the local kappa_d.
            {
                auto const& ehlers_material =
                    static_cast<MaterialLib::Solids::Ehlers::SolidEhlers<
                        DisplacementDim> const&>(_ip_data[ip].solid_material);
                auto const damage_properties =
                    ehlers_material.evaluatedDamageProperties(t, x_position);
                auto const material_properties =
                    ehlers_material.evaluatedMaterialProperties(t, x_position);

                // Ehlers material state variables
                auto& state_vars =
                    static_cast<MaterialLib::Solids::Ehlers::StateVariables<
                        DisplacementDim>&>(
                        *_ip_data[ip].material_state_variables);

                double const eps_p_eff_diff =
                    state_vars.eps_p.eff - state_vars.eps_p_prev.eff;

                _ip_data[ip].kappa_d = calculateDamageKappaD<DisplacementDim>(
                    eps_p_eff_diff, sigma, _ip_data[ip].kappa_d_prev,
                    damage_properties.h_d, material_properties);

                if (!_ip_data[ip].active_self)
                {
                    _ip_data[ip].active_self |= _ip_data[ip].kappa_d > 0;
                    if (_ip_data[ip].active_self)
                    {
                        for (auto const& tuple :
                             _ip_data[ip].non_local_assemblers)
                        {
                            // Activate the integration point.
                            tuple.ip_l_pointer->activated = true;
                        }
                    }
                }
            }
        }
    }

    void assembleWithJacobian(double const t, double const /*dt*/,
                              std::vector<double> const& local_x,
                              std::vector<double> const& /*local_xdot*/,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_b_data,
                              std::vector<double>& local_Jac_data) override
    {
        auto const local_matrix_size = local_x.size();

        auto local_Jac = MathLib::createZeroedMatrix<StiffnessMatrixType>(
            local_Jac_data, local_matrix_size, local_matrix_size);

        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            local_b_data, local_matrix_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition x_position;
        x_position.setElementID(_element.getID());

        // Non-local integration.
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            x_position.setIntegrationPoint(ip);
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(_element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);

            auto& sigma = _ip_data[ip].sigma;
            auto& C = _ip_data[ip].C;
            double& damage = _ip_data[ip].damage;

            {
                double nonlocal_kappa_d = 0;

                if (_ip_data[ip].active_self || _ip_data[ip].activated)
                {
                    for (auto const& tuple : _ip_data[ip].non_local_assemblers)
                    {
                        // Get local variable for the integration point l.
                        double const kappa_d_l = tuple.ip_l_pointer->kappa_d;
                        double const a_kl_times_w_l = tuple.alpha_kl_times_w_l;
                        nonlocal_kappa_d += a_kl_times_w_l * kappa_d_l;
                    }
                }

                auto const& ehlers_material =
                    static_cast<MaterialLib::Solids::Ehlers::SolidEhlers<
                        DisplacementDim> const&>(_ip_data[ip].solid_material);

                //
                // Overnonlocal formulation
                //
                // See (Di Luzio & Bazant 2005, IJSS) for details.
                // The implementation would go here and would be for a given
                // gamma_nonlocal:
                //
                // Update nonlocal damage with local damage (scaled with 1 -
                // \gamma_{nonlocal}) for the current integration point and the
                // nonlocal integral part.
                // nonlocal_kappa_d = (1. - gamma_nonlocal) * kappa_d +
                //                    gamma_nonlocal * nonlocal_kappa_d;

                nonlocal_kappa_d = std::max(0., nonlocal_kappa_d);

                // Update damage based on nonlocal kappa_d
                {
                    auto const damage_properties =
                        ehlers_material.evaluatedDamageProperties(t,
                                                                  x_position);
                    damage = calculateDamage(nonlocal_kappa_d,
                                             damage_properties.alpha_d,
                                             damage_properties.beta_d);
                    damage = std::max(0., damage);
                }
                sigma = sigma * (1. - damage);
            }

            local_b.noalias() -= B.transpose() * sigma * w;
            local_Jac.noalias() += B.transpose() * C * (1. - damage) * B * w;
        }
    }

    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void postTimestepConcrete(Eigen::VectorXd const& /*local_x*/,
                              double const /*t*/, double const /*dt*/) override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data[ip].pushBackState();
        }
    }

    void computeCrackIntegral(std::size_t mesh_item_id,
                              NumLib::LocalToGlobalIndexMap const& dof_table,
                              GlobalVector const& x,
                              double& crack_volume) override
    {
        auto const indices = NumLib::getIndices(mesh_item_id, dof_table);
        auto local_x = x.get(indices);

        auto u = Eigen::Map<typename BMatricesType::NodalForceVectorType const>(
            local_x.data(), ShapeFunction::NPOINTS * DisplacementDim);

        int const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (int ip = 0; ip < n_integration_points; ip++)
        {
            auto const& dNdx = _ip_data[ip].dNdx;
            auto const& d = _ip_data[ip].damage;
            auto const& w = _ip_data[ip].integration_weight;

            double const div_u =
                Deformation::divergence<DisplacementDim,
                                        ShapeFunction::NPOINTS>(u, dNdx);
            crack_volume += div_u * d * w;
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _secondary_data.N[integration_point];

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    std::vector<double> const& getNodalValues(
        std::vector<double>& nodal_values) const override
    {
        nodal_values.clear();
        auto local_b = MathLib::createZeroedVector<NodalDisplacementVectorType>(
            nodal_values, ShapeFunction::NPOINTS * DisplacementDim);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& w = _ip_data[ip].integration_weight;

            auto const& N = _ip_data[ip].N;
            auto const& dNdx = _ip_data[ip].dNdx;

            auto const x_coord =
                NumLib::interpolateXCoordinate<ShapeFunction,
                                               ShapeMatricesType>(_element, N);
            auto const B = LinearBMatrix::computeBMatrix<
                DisplacementDim, ShapeFunction::NPOINTS,
                typename BMatricesType::BMatrixType>(dNdx, N, x_coord,
                                                     _is_axially_symmetric);
            auto& sigma = _ip_data[ip].sigma;

            local_b.noalias() += B.transpose() * sigma * w;
        }

        return nodal_values;
    }

    std::vector<double> const& getIntPtFreeEnergyDensity(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        transform(
            cbegin(_ip_data), cend(_ip_data), back_inserter(cache),
            [](auto const& ip_data) { return ip_data.free_energy_density; });

        return cache;
    }

    std::vector<double> const& getIntPtEpsPV(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        transform(cbegin(_ip_data), cend(_ip_data), back_inserter(cache),
                  [](auto const& ip_data) { return *ip_data.eps_p_V; });

        return cache;
    }

    std::vector<double> const& getIntPtEpsPDXX(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        transform(cbegin(_ip_data), cend(_ip_data), back_inserter(cache),
                  [](auto const& ip_data) { return *ip_data.eps_p_D_xx; });

        return cache;
    }

    std::vector<double> const& getIntPtSigma(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma, cache);
    }

    std::vector<double> const& getIntPtEpsilon(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::eps, cache);
    }

    std::size_t setSigma(double const* values)
    {
        return ProcessLib::setIntegrationPointKelvinVectorData<DisplacementDim>(
            values, _ip_data, &IpData::sigma);
    }

    // TODO (naumov) This method is same as getIntPtSigma but for arguments and
    // the ordering of the cache_mat.
    // There should be only one.
    std::vector<double> getSigma() const override
    {
        return ProcessLib::getIntegrationPointKelvinVectorData<DisplacementDim>(
            _ip_data, &IpData::sigma);
    }

    void setKappaD(double value)
    {
        for (auto& ip_data : _ip_data)
        {
            ip_data.kappa_d = value;
        }
    }
    std::vector<double> getKappaD() const override
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        std::vector<double> result_values;
        result_values.resize(n_integration_points);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            result_values[ip] = _ip_data[ip].kappa_d;
        }

        return result_values;
    }

    std::vector<double> const& getIntPtDamage(
        const double /*t*/,
        std::vector<GlobalVector*> const& /*x*/,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& /*dof_table*/,
        std::vector<double>& cache) const override
    {
        return ProcessLib::getIntegrationPointScalarData(
            _ip_data, &IpData::damage, cache);
    }

    unsigned getNumberOfIntegrationPoints() const override
    {
        return _integration_method.getNumberOfPoints();
    }

    typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(int const integration_point) const override
    {
        return *_ip_data[integration_point].material_state_variables;
    }

private:
    std::vector<double> const& getIntPtSigma(std::vector<double>& cache,
                                             std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)
            {  // xx, yy, zz components
                cache.push_back(ip_data.sigma[component]);
            }
            else
            {  // mixed xy, yz, xz components
                cache.push_back(ip_data.sigma[component] / std::sqrt(2));
            }
        }

        return cache;
    }

    std::vector<double> const& getIntPtEpsilon(
        std::vector<double>& cache, std::size_t const component) const
    {
        cache.clear();
        cache.reserve(_ip_data.size());

        for (auto const& ip_data : _ip_data)
        {
            if (component < 3)  // xx, yy, zz components
                cache.push_back(ip_data.eps[component]);
            else  // mixed xy, yz, xz components
                cache.push_back(ip_data.eps[component] / std::sqrt(2));
        }

        return cache;
    }

    IntegrationPointDataNonlocalInterface*
    getIPDataPtr(int const ip) override
    {
        return &_ip_data[ip];
    }

private:
    SmallDeformationNonlocalProcessData<DisplacementDim>& _process_data;

    std::vector<IpData, Eigen::aligned_allocator<IpData>> _ip_data;

    IntegrationMethod const _integration_method;
    MeshLib::Element const& _element;
    SecondaryData<typename ShapeMatrices::ShapeType> _secondary_data;
    bool const _is_axially_symmetric;

    static const int displacement_size =
        ShapeFunction::NPOINTS * DisplacementDim;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
