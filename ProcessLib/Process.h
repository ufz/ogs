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
    /// for MPI based parallel computing.
    /// Since the global matrix, vector and linear equation in ProjectData
    /// are defined as smarter point type (unique_ptr or might be shared_ptr)
    /// variables, their memory occupations are automatically released at
    /// the end of ProjectData terminated, i.e. the end of the main program.
    /// However if  the global matrix, vector and linear equation are created
    /// under MPI environment, their memory occupations must be released before
    /// calling of MPI_Finalize, PetscFinalize. That is why the following fucntion
    /// is needed.
    virtual void releaseEquationMemory() = 0;

    /// Postprocessing after solve().
    /// The file_name is indicating the name of possible output file.
    virtual void post(std::string const& file_name) = 0;
    virtual void postTimestep(std::string const& file_name, const unsigned timestep) = 0;

protected:
    MeshLib::Mesh& _mesh;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_H_
