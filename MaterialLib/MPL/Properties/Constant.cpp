/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 7, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Constant.h"

namespace MaterialPropertyLib
{
struct ZeroInitPropertyDataType
{
    PropertyDataType operator()(double) const { return 0.; }

    PropertyDataType operator()(Eigen::Vector2d) const
    {
        return Eigen::Vector2d::Zero().eval();
    }

    PropertyDataType operator()(Eigen::Vector3d) const
    {
        return Eigen::Vector3d::Zero().eval();
    }

    PropertyDataType operator()(Eigen::Matrix<double, 2, 2>) const
    {
        return Eigen::Matrix<double, 2, 2>::Zero().eval();
    }
    PropertyDataType operator()(Eigen::Matrix<double, 3, 3>) const
    {
        return Eigen::Matrix<double, 3, 3>::Zero().eval();
    }

    PropertyDataType operator()(Eigen::Matrix<double, 4, 1>) const
    {
        return Eigen::Matrix<double, 4, 1>::Zero().eval();
    }

    PropertyDataType operator()(Eigen::Matrix<double, 6, 1>) const
    {
        return Eigen::Matrix<double, 6, 1>::Zero().eval();
    }

    PropertyDataType operator()(Eigen::MatrixXd) const
    {
        return Eigen::MatrixXd(0, 0);
    }
};

Constant::Constant(std::string name, PropertyDataType const& v)
{
    name_ = std::move(name);
    value_ = v;
    dvalue_ = std::visit(ZeroInitPropertyDataType{}, v);
};
}  // namespace MaterialPropertyLib
