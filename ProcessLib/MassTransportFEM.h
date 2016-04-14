/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_MASSTRANSPORT_FEM_H_
#define PROCESS_LIB_MASSTRANSPORT_FEM_H_

#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "Parameter.h"
#include "ProcessUtil.h"


namespace ProcessLib
{

namespace MassTransport
{

// TODO now this interface is basically the same for all processes that assemble a
//      FirstOrderImplicitQuasiLinear ODE system.
template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&,
            GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    LocalAssemblerData(MeshLib::Element const& e,
                       std::size_t const local_matrix_size,
                       unsigned const integration_order,
                       Parameter<double, MeshLib::Element const&> const&
						diffusion_coefficient,
						Parameter<double, MeshLib::Element const&> const&
						velocity)//std::vector<double>
        : _shape_matrices(
              initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
                  e, integration_order))
        , _diffusion_coefficient([&diffusion_coefficient, &e]()
          {
              return diffusion_coefficient(e);
          })
		, _velocity([&velocity, &e]()
		  {
			  return velocity(e);
		  })
        // TODO narrowing conversion
		, _localM(local_matrix_size, local_matrix_size)
        , _localA(local_matrix_size, local_matrix_size)
        , _localRhs(local_matrix_size)
        , _integration_order(integration_order)
    {}

    void assemble(double const /*t*/, std::vector<double> const& /*local_x*/) override
    {
        _localA.setZero();
		_localM.setZero();
        _localRhs.setZero();
		double rho = 1.;
		//double test = _velocity();
		Eigen::VectorXd vel = Eigen::VectorXd::Zero(GlobalDim);
		for (int ii = 0; ii < GlobalDim; ii++) {
			vel(ii) = _velocity();
		}
        IntegrationMethod_ integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();//retuen gauss point number

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            _localA.noalias() += sm.dNdx.transpose() *
                                 _diffusion_coefficient() * sm.dNdx *
                                 sm.detJ * wp.getWeight();
			/*_localA.noalias() += sm.N.transpose() *
				Eigen::Map<Eigen::Vector3d>(_velocity().data()) * sm.dNdx *
				sm.detJ * wp.getWeight();*/
			_localA.noalias() += sm.N *vel.transpose()* sm.dNdx *
				sm.detJ * wp.getWeight();
			_localM.noalias() += sm.N *
				rho * sm.N.transpose() *
				sm.detJ * wp.getWeight();//Eigen::Map<Eigen::VectorXd>(
        }
    }

    void addToGlobal(AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
        const override
    {
		M.add(indices, _localM);
        K.add(indices, _localA);
        b.add(indices.rows, _localRhs);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    std::function<double(void)> _diffusion_coefficient;
	std::function<double(void)> _velocity;
	//std::function<std::vector<double>(void)> _velocity;
    NodalMatrixType _localA;
	NodalMatrixType _localM;
	//NodalMatrixType _localD;
    NodalVectorType _localRhs;

    unsigned const _integration_order;
};


}   // namespace MassTransport
}   // namespace ProcessLib

#endif  // PROCESS_LIB_MASSTRANSPORT_FEM_H_
