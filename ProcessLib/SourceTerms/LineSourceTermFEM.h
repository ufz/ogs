/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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

const unsigned NUM_NODAL_DOF = 1;

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LineSourceTermLocalAssembler final
    : public LineSourceTermLocalAssemblerInterface
{
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
        : _parameter(line_source_term_parameter),
          _integration_method(integration_order),
          _local_rhs(local_matrix_size)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                _integration_method.getWeightedPoint(ip).getWeight() *
                shape_matrices[ip].integralMeasure * shape_matrices[ip].detJ *
                shape_matrices[ip].N);
        }
    }

    void integrate(std::size_t const id,
                   NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
                   double const t, GlobalVector& b) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(id);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            pos.setIntegrationPoint(ip);
            auto const st_val = _parameter(t, pos)[0];

            _local_rhs.noalias() += st_val * _ip_data[ip];
        }
        auto const indices = NumLib::getIndices(id, source_term_dof_table);
        b.add(indices, _local_rhs);
    }

private:
    ParameterLib::Parameter<double> const& _parameter;

    IntegrationMethod const _integration_method;
    std::vector<NodalRowVectorType,
                Eigen::aligned_allocator<NodalRowVectorType>>
        _ip_data;
    NodalVectorType _local_rhs;
};

}  // namespace ProcessLib
