#include <cassert>
#include <cstdio>

#include "logog/include/logog.hpp"

#include "TESProcess-notpl.h"
#include "TESFEM-data-fwd.h"



namespace ProcessLib
{

namespace TES
{

void printGlobalMatrix(const Eigen::SparseMatrix<double> &mat)
{
    const unsigned R = mat.rows();
    const unsigned C = mat.cols();
    assert(R==C);


#if 0
    const unsigned N = R/NODAL_DOF;

    for (unsigned nr=0; nr<N; ++nr)
    {
        for (unsigned vr=0; vr<NODAL_DOF; ++vr)
        {
            const unsigned row = vr*N+nr;

            for (unsigned nc=0; nc<N; ++nc)
            {
                for (unsigned vc=0; vc<NODAL_DOF; ++vc)
                {
                    const unsigned col = vc*N+nc;
                    if (col!=0) std::printf("  ");

                    std::printf("%14g", const_cast<Eigen::SparseMatrix<double>& >(mat).coeffRef(row, col));
                    // std::cout << r << ", " << c << ": " << _A->getRawMatrix().coeffRef(r,c) << std::endl;
                }
            }

            std::printf("\n");
        }
    }
#else
    for (unsigned r=0; r<R; ++r)
    {
        for (unsigned c=0; c<C; ++c)
        {
            std::printf("%19.12g", const_cast<Eigen::SparseMatrix<double>& >(mat).coeffRef(r, c));
        }
        std::printf("\n");
    }
#endif
}

void printGlobalVector(const Eigen::Ref<Eigen::VectorXd>& vec)
{
    for (unsigned i=0; i<vec.size(); ++i)
    {
        std::printf("%19.12g\n", vec[i]);
    }
}



bool enforceConstraint(Eigen::VectorXd* vec,
                       const Eigen::Ref<Eigen::VectorXd>& previous_solution,
                       const unsigned dof, const unsigned numDof, const unsigned nnodes,
                       const double vmin, const double vmax)
{
    assert(vec->size() == nnodes*numDof);

    bool ok = true;

    for (std::size_t i=0; i<nnodes; ++i)
    {
        auto& v = (*vec)[dof*nnodes+i]; // TODO [CL] depends on local<->global mapping
        if (vmin > v) {
            v = 0.5 * (previous_solution[dof*nnodes+i] + vmin);
            ok = false;
        } else if (vmax < v) {
            v = 0.5 * (previous_solution[dof*nnodes+i] + vmax);
            ok = false;
        }
    }

    return ok;
}



bool calculateError(Eigen::VectorXd* current_solution,
                    const Eigen::Ref<Eigen::VectorXd>& previous_solution,
                    AssemblyParams* materials)
{
    auto const& current_s = *current_solution;
    assert(current_s.size() == previous_solution.size());

    double rel_err = 0.0;
    auto const N = current_s.size();
    auto const dmin = std::numeric_limits<double>::min();
    auto const deps = std::numeric_limits<double>::epsilon();

    for (std::remove_const<decltype(N)>::type i=0; i<N; ++i)
    {
        // const double c = const_cast<Eigen::SparseMatrix<double>&>(current_solution).coeffRef(i,1);
        // const double p = const_cast<Eigen::SparseMatrix<double>&>(previous_solution).coeffRef(i,1);
        const double c = current_s[i];
        const double p = previous_solution[i];
        if (std::fabs(p) < 10.0 * dmin)
        { // previous value is almost equal to zero
            if (std::fabs(c-p) >= 10.0 * deps)
            { // values are not almost equal
                // TODO [CL] what about very large c?
                rel_err += std::fabs(c-p);
            }
            else
            { // values are almost equal -- error of zero
            }
        }
        else
        {
            // previous value is different from zero

            rel_err += std::fabs((c - p)/p);
        }
    }

    rel_err /= N;

    DBUG("average relative error per d.o.f.: %e", rel_err);


    std::puts("-------- solution -------\n");
    printGlobalVector(* current_solution);

    // check constraints on the solution

    const bool constr_ok0 = enforceConstraint(current_solution, previous_solution, 0, NODAL_DOF, N/NODAL_DOF, 0.0, std::numeric_limits<double>::max()); // pressure
    const bool constr_ok1 = enforceConstraint(current_solution, previous_solution, 1, NODAL_DOF, N/NODAL_DOF, 0.0, std::numeric_limits<double>::max()); // temperature
    const bool constr_ok2 = enforceConstraint(current_solution, previous_solution, 2, NODAL_DOF, N/NODAL_DOF, 0.0, 1.0); // concentration

    bool constr_ok = constr_ok0 && constr_ok1 && constr_ok2;

    if (! constr_ok)
    {
        // decrease time step
        const double old_ts = materials->_delta_t;
        const double new_ts = old_ts / 2.0;
        DBUG("some constraints were violated. reducing timestep from %g to %g", old_ts, new_ts);
        // materials->_time_step = new_ts;

        return false;
    }

    return rel_err < NONLINEAR_ERROR_TOLERANCE;
}

} // namespace TES

} // namespace ProcessLib
