/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <set>
#include <string>
#include <vector>

#include "ProcessLib/Output/OutputDataSpecification.h"

namespace ProcessLib
{
enum class OutputType : uint8_t
{
    vtk,
    xdmf
};

struct OutputConfig
{
    /// Supported output types are "VTK" and "XDMF".
    ///
    /// For the "VTK" type a `.pvd` file referencing each time step, which is
    /// stored in a `.vtu` file is written. For the "XDMF" type all time steps
    /// are written to a `.h5` HDF5 file and a `.xdmf` file references the
    /// actual data in the `.h5` file.  Both, the `.pvd` and the `.xdmf` files
    /// can be loaded by Paraview and other VTK based post-processing tools.
    OutputType output_type;
    std::string prefix;
    std::string suffix;
    bool compress_output;
    unsigned int number_of_files;
    unsigned int chunk_size_bytes;
    std::string data_mode;
    /// A list of repeat/step-count pairs. If the list is empty, and no
    /// fixed_output_times were specified, a default pair 1/1 will be inserted
    /// resulting in all time steps being written.
    ///
    /// The last pair is repeated until end of simulation time is reached.
    std::vector<PairRepeatEachSteps> repeats_each_steps;
    /// A list of variable names for output. If the list is empty, all available
    /// variables are written.
    ///
    /// The variable names include: primary variables like `temperature`,
    /// `displacement`, *etc.*, secondary variables like `epsilon` or `sigma`,
    /// residuals `HeatFlux`, `NodalForces`, integration point data (required
    /// for restart) for example `sigma_ip`.  Available names are process and
    /// material models specific and can not be listed here.
    std::set<std::string> output_variables;
    bool output_extrapolation_residuals;
    std::vector<std::string> mesh_names_for_output;
    /// A list of points in time for output. These fixed output times are taken
    /// by the simulation independent of what time stepping scheme is specified.
    ///
    /// Often used to arrive at some critical points in the simulation like load
    /// change.
    std::vector<double> fixed_output_times;
    bool output_iteration_results;
};
}  // namespace ProcessLib
