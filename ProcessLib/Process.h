/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include <string>

namespace MeshLib
{
    class Mesh;
}

namespace ProcessLib
{

class Process
{
public:
    Process(MeshLib::Mesh& mesh)
        : _mesh(mesh)
    { }

    virtual ~Process() = default;

    virtual void initialize() = 0;
    virtual bool solve(const double delta_t) = 0;

    /// Explicitly release memory of global vector, matrix and linear solvers
    /// for MPI based parallel computing
    virtual void releaseEquationMemory()
    {
        // This is dummy code to fix _compilation_ of currently broken petsc
        // configuration.
    }

    /// Postprocessing after solve().
    /// The file_name is indicating the name of possible output file.
    virtual void post(std::string const& file_name) = 0;
    virtual void postTimestep(std::string const& file_name, const unsigned timestep) = 0;

protected:
    MeshLib::Mesh& _mesh;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
