/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Mesh.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
{
struct FlagAndFluxAndDFlux
{
    //! Indicates if flux and dflux have been computed.
    //! If false, don't use FlagAndFluxAndDFlux#flux nor
    //! FlagAndFluxAndDFlux#dFlux.
    bool flag;

    //! The computed flux.
    double flux;

    //! The partial derivates of FlagAndFluxAndDFlux#flux w.r.t. each primary
    //! variable.
    std::vector<double> dFlux;
};

//! Contains data commonly used by Python BCs and STs, in particular by their
//! local assemblers.
template <typename BcOrStPythonSideInterface>
struct BcOrStData
{
    //! Python object computing BC or ST values.
    BcOrStPythonSideInterface const* const bc_or_st_object;

    //! Global component ID of the (variable, component) to which this BC or ST
    //! is applied.
    int const global_component_id;

    //! The domain, where this BC or ST will be applied.
    MeshLib::Mesh const& bc_or_st_mesh;

    //! This BC or ST is applied to the global equation system of some process.
    //! These are all process variables of that process.
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process;

    //! The shape function order of the process variable to which this BC or ST
    //! belongs.
    unsigned const shape_function_order;
};

}  // namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
