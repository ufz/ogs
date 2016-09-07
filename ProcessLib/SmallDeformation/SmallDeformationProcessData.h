/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_SMALLDEFORMATIONPROCESSDATA_H_
#define PROCESSLIB_SMALLDEFORMATION_SMALLDEFORMATIONPROCESSDATA_H_

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
struct SmallDeformationProcessData
{
    SmallDeformationProcessData(
        std::unique_ptr<Solids::MechanicsBase<DisplacementDim>>&& material)
        : _material{std::move(material)}
    {
    }

    SmallDeformationProcessData(SmallDeformationProcessData&& other)
        : _material{std::move(other._material)}
    {
    }

    //! Copies are forbidden.
    SmallDeformationProcessData(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(SmallDeformationProcessData&&) = delete;

    std::unique_ptr<Solids::MechanicsBase<DisplacementDim>> _material;
    double dt;
    double t;
};

}  // namespace SmallDeformation
}  // namespace ProcessLib

#endif  // PROCESSLIB_SMALLDEFORMATION_SMALLDEFORMATIONPROCESSDATA_H_
