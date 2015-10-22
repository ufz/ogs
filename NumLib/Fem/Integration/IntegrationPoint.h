/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_INTEGRATIONPOINT_H_
#define NUMLIB_INTEGRATIONPOINT_H_

namespace NumLib
{
/// Integration rule for point elements.
///
/// The integration order is not stored or used for point integration. It is
/// only needed to satisfy the common integration rule concepts.
class IntegrationPoint
{
	typedef MathLib::TemplateWeightedPoint<double, double, 1> WeightedPoint;

public:
	/// IntegrationPoint constructor for given order.
	explicit IntegrationPoint(std::size_t /* order */)
	{
	}

	/// Change the integration order.
	void setIntegrationOrder(std::size_t /* order */)
	{
	}

	/// return current integration order.
	std::size_t getIntegrationOrder() const
	{
		return 0;
	}

	/// return the number of sampling points
	std::size_t getNPoints() const
	{
		return 1;
	}

	/**
	 * get coordinates of a integration point
	 *
	 * @param igp      The integration point index
	 * @return a weighted point
	 */
	WeightedPoint getWeightedPoint(std::size_t igp)
	{
		return getWeightedPoint(getIntegrationOrder(), igp);
	}

	/**
	 * get coordinates of a integration point
	 *
	 * @param order    the number of integration points
	 * @param pt_id     the sampling point id
	 * @return weight
	 */
	static WeightedPoint getWeightedPoint(std::size_t /* order */, std::size_t /*igp*/)
	{
		return WeightedPoint({{1}}, 1);
	}

	template <typename Method>
	static WeightedPoint getWeightedPoint(std::size_t /*igp*/)
	{
		return WeightedPoint({{1}}, 1);
	}

	/**
	 * get the number of integration points
	 *
	 * @param order    the number of integration points
	 * @return the number of points
	 */
	static std::size_t getNPoints(std::size_t /* order */)
	{
		return 1;
	}
};
}

#endif  // NUMLIB_INTEGRATIONPOINT_H_
