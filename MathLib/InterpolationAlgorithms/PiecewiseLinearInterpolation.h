/**
 * \file PiecewiseLinearInterpolation.h
 *
 * Created on 2010-09-07 by Thomas Fischer
 */

#ifndef PIECEWISELINEARINTERPOLATION_H_
#define PIECEWISELINEARINTERPOLATION_H_

#include <limits>
#include <vector>

namespace MathLib
{
class PiecewiseLinearInterpolation
{
public:
	PiecewiseLinearInterpolation(const std::vector<double>& supporting_points,
	                    const std::vector<double>& values_at_supp_pnts);
	PiecewiseLinearInterpolation(const std::vector<double>& supporting_points,
	                    const std::vector<double>& values_at_supp_pnts,
	                    const std::vector<double>& points_to_interpolate,
	                    std::vector<double>& values_at_interpol_pnts);
	virtual ~PiecewiseLinearInterpolation();

	double getValue ( double pnt_to_interpolate );

private:
	const std::vector<double>& _supporting_points;
	const std::vector<double>& _values_at_supp_pnts;
};
} // end namespace MathLib

#endif /* PIECEWISELINEARINTERPOLATION_H_ */
