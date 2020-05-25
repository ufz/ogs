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

#include <vector>

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
class LineSourceTermLocalAssemblerInterface
{
public:
    virtual void integrate(
        std::size_t const id,
        NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
        double const t, GlobalVector& b) = 0;
    virtual ~LineSourceTermLocalAssemblerInterface() = default;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LineSourceTermLocalAssembler final
    : public LineSourceTermLocalAssemblerInterface
{
    static const unsigned NUM_NODAL_DOF = 1;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

    using LocalAssemblerTraits = ProcessLib::LocalAssemblerTraits<
        ShapeMatricesType, ShapeFunction::NPOINTS, NUM_NODAL_DOF, GlobalDim>;

    using NodalVectorType = typename LocalAssemblerTraits::LocalVector;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

public:
    LineSourceTermLocalAssembler(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        ParameterLib::Parameter<double> const& line_source_term_parameter)
        : parameter_(line_source_term_parameter),
          integration_method_(integration_order),
          local_rhs_(local_matrix_size)
    {
        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, integration_method_);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            ip_data_.emplace_back(
                integration_method_.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ *
                shape_matrices[ip].N);
        }
    }

    void integrate(std::size_t const id,
                   NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
                   double const t, GlobalVector& b) override
    {
        local_rhs_.setZero();

        unsigned const n_integration_points =
            integration_method_.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(id);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const st_val = parameter_(t, pos)[0];

            local_rhs_.noalias() += st_val * ip_data_[ip];
        }
        auto const indices = NumLib::getIndices(id, source_term_dof_table);
        b.add(indices, local_rhs_);
    }

private:
    ParameterLib::Parameter<double> const& parameter_;

    IntegrationMethod const integration_method_;
    std::vector<NodalRowVectorType,
                Eigen::aligned_allocator<NodalRowVectorType>>
        ip_data_;
    NodalVectorType local_rhs_;
};

}  // namespace ProcessLib
