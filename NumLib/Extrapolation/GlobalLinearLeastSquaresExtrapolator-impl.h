#include "logog/include/logog.hpp"

#include "Eigen/SparseQR"

#include "Util.h"

#include "GlobalLinearLeastSquaresExtrapolator.h"

namespace
{

template<typename GlobalMatrix, typename GlobalVector,
         typename VariableEnum, typename LocalAssembler>
void
GLLSQ_gatherElementData(
        LocalAssembler const* loc_asm, VariableEnum var,
        AssemblerLib::LocalToGlobalIndexMap::LineIndex const& indices,
        GlobalMatrix& matrix, GlobalVector& integration_point_values,
        std::size_t& start_row
        )
{
    auto const& gp_vals = loc_asm->getIntegrationPointValues(var);

    const unsigned nn = loc_asm->getShapeMatrix(0).rows(); // number of mesh nodes
    const unsigned ni = gp_vals.size();                    // number of integration points

    for (unsigned gp=0; gp<ni; ++gp)
    {
        integration_point_values[start_row + gp] = gp_vals[gp];

        auto const& shp_mat = loc_asm->getShapeMatrix(gp);
        assert(shp_mat.rows() == nn);

        for (unsigned n=0; n<nn; ++n)
        {
            matrix.setValue(start_row+gp, indices[n], shp_mat(n));
        }
    }

    start_row += ni;
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
double
GLLSQ_calculateResiudalElement(LocalAssembler const* loc_asm, VariableEnum var,
        AssemblerLib::LocalToGlobalIndexMap::LineIndex const& indices,
        GlobalVector const& nodal_vals)
{
    auto const& gp_vals = loc_asm->getIntegrationPointValues(var);
    const unsigned ni = gp_vals.size();        // number of gauss points

    // filter nodal values of the current element
    std::vector<double> nodal_vals_element;
    nodal_vals_element.resize(indices.size());
    for (unsigned i=0; i<indices.size(); ++i) {
        nodal_vals_element[i] = nodal_vals[indices[i]];
    }

    double residual = 0.0;
    double gp_val_extrapol = 0.0;
    double* gp_val_extrapol2[1] = { &gp_val_extrapol };

    for (unsigned gp=0; gp<ni; ++gp)
    {
        ProcessLib::shapeFunctionInterpolate(
                    nodal_vals_element, loc_asm->getShapeMatrix(gp),
                    1, gp_val_extrapol2);
        auto const& ax_m_b = gp_val_extrapol - gp_vals[gp];
        residual += ax_m_b * ax_m_b;
    }

    return residual / ni;
}


template<typename Matrix, typename RHS, typename Unknowns>
void
GLLSQ_solveLinearLeastSquares(Matrix const& mat, RHS const& rhs, Unknowns& unknowns)
{
    switch(ProcessLib::linear_least_squares)
    {
    case ProcessLib::LinearLeastSquaresBy::NormalEquation:
    {
        DBUG("solving normal equation...");
        // unknowns = (mat.transpose() * mat).ldlt().solve(mat.transpose() * rhs);
        break;
    }
    case ProcessLib::LinearLeastSquaresBy::QR:
    {
        Eigen::SparseQR<Matrix, Eigen::COLAMDOrdering<int>> solver;
        solver.compute(mat);
        unknowns = solver.solve(rhs);
    }
    default:
        ERR("chosen linear least squares method not yet implemented.");
    }
}

} // anonymous namespace


namespace ProcessLib
{

template<typename GlobalMatrix, typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
GlobalLinearLeastSquaresExtrapolator<GlobalMatrix, GlobalVector, VariableEnum, LocalAssembler>::
extrapolate(LocalAssemblers const& loc_asms, VariableEnum var)
{
    _nodal_values = 0.0;

    std::size_t num_integration_points = 0;
    std::size_t num_integration_points_per_element = 0;

    for (std::size_t ei=0; ei<loc_asms.size(); ++ei)
    {
        auto const n = loc_asms[ei]->getIntegrationPointValues(var).size();
        num_integration_points += n;
        if (n > num_integration_points_per_element) num_integration_points_per_element = n;
    }

    GlobalMatrix mat(num_integration_points, num_integration_points_per_element);
    _integration_point_values.resize(num_integration_points);

    std::size_t start_row = 0;

    for (std::size_t i=0; i<loc_asms.size(); ++i)
    {
        GLLSQ_gatherElementData(loc_asms[i], var, _local_to_global[i].rows,
                          mat, _integration_point_values, start_row);
    }

    GLLSQ_solveLinearLeastSquares(mat.getRawMatrix(),
                                  _integration_point_values.getRawVector(),
                                  _nodal_values.getRawVector());
}

template<typename GlobalMatrix, typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
GlobalLinearLeastSquaresExtrapolator<GlobalMatrix, GlobalVector, VariableEnum, LocalAssembler>::
calculateResiduals(LocalAssemblers const& loc_asms, VariableEnum var)
{
    _residuals.resize(loc_asms.size());

    for (std::size_t ei=0; ei<loc_asms.size(); ++ei)
    {
        _residuals[ei] = GLLSQ_calculateResiudalElement(
                    loc_asms[ei], var, _local_to_global[ei].rows,
                    _nodal_values
                    );
    }
}

} // namespace ProcessLib
