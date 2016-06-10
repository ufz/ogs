/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_LOCALASSEMBLERINTERFACE_H
#define PROCESSLIB_LOCALASSEMBLERINTERFACE_H

namespace ProcessLib
{

/*! Common interface for local assemblers
 * NumLib::ODESystemTag::FirstOrderImplicitQuasilinear ODE systems.
 *
 * \todo Generalize to other NumLib::ODESystemTag's.
 */
template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerInterface
{
public:
    virtual ~LocalAssemblerInterface() = default;

    virtual void assemble(double const t, std::vector<double> const& local_x) = 0;

    virtual void addToGlobal(
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const&,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const = 0;

    virtual void assembleJacobian(double const /*t*/,
                                  std::vector<double> const& /*local_x*/)
    {
        OGS_FATAL(
            "assembleJacobian function is not implemented in the local "
            "assembler.");
    }

    virtual void addJacobianToGlobal(NumLib::LocalToGlobalIndexMap::
                                         RowColumnIndices const& /*indices*/,
                                     GlobalMatrix& /*Jac*/) const
    {
        OGS_FATAL(
            "addJacobianToGlobal function is not implemented in the local "
            "assembler.");
    }

    virtual void preTimestep(std::vector<double> const& /*local_x*/,
                             double const /*t*/, double const /*delta_t*/)
    {
    }
    virtual void postTimestep(std::vector<double> const& /*local_x*/) {}
};

} // namespace ProcessLib

#endif // PROCESSLIB_LOCALASSEMBLERINTERFACE_H
