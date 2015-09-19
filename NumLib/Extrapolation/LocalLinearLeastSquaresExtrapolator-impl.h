#include "logog/include/logog.hpp"

#include "Util.h"

#include "LocalLinearLeastSquaresExtrapolator.h"


// see http://eigen.tuxfamily.org/dox-devel/group__LeastSquares.html
enum class LinearLeastSquaresBy { SVD, QR, NormalEquation };
const LinearLeastSquaresBy linear_least_squares = LinearLeastSquaresBy::NormalEquation;


namespace
{
template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LLLSQ_extrapolateElement(
        LocalAssembler const* loc_asm, VariableEnum var,
        AssemblerLib::LocalToGlobalIndexMap::LineIndex const& indices,
        GlobalVector& nodal_vals, GlobalVector& counts
        )
{
    // auto const& shp_mats = loc_asm->getShapeMatrices();
    auto const& gp_vals = loc_asm->getIntegrationPointValues(var);

    const unsigned nn = loc_asm->getShapeMatrix(0).rows(); // number of mesh nodes
    // const unsigned nn = 5;
    const unsigned ni = gp_vals.size();        // number of gauss points


    Eigen::MatrixXd N(ni, nn);

    for (unsigned gp=0; gp<ni; ++gp)
    {
        auto const& shp_mat = loc_asm->getShapeMatrix(gp);
        assert(shp_mat.rows() == nn);

        for (unsigned n=0; n<nn; ++n)
        {
            N(gp, n) = shp_mat(n);
        }
    }

    const Eigen::Map<const Eigen::VectorXd> gpvs(gp_vals.data(), gp_vals.size());

    switch (linear_least_squares)
    {
    case LinearLeastSquaresBy::NormalEquation:
    {
        DBUG("solving normal equation...");
        Eigen::VectorXd elem_nodal_vals = (N.transpose() * N).ldlt().solve(N.transpose() * gpvs);
        nodal_vals.add(indices, elem_nodal_vals);
        counts.add(indices, 1.0);
        break;
        // return (N.transpose() * N).ldlt().solve(N.transpose() * gpvs);
    }
    default:
        ERR("chosen linear least squares method not yet implemented.");
    }
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
double
calculateResiudalElement(LocalAssembler const* loc_asm, VariableEnum var,
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

} // anonymous namespace


namespace ProcessLib
{

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
extrapolate(LocalAssemblers const& loc_asms, VariableEnum var)
{
    _nodal_values = 0.0;
    // _counts = 0.0;

    GlobalVector counts(_nodal_values.size());

    for (std::size_t i=0; i<loc_asms.size(); ++i)
    {
        LLLSQ_extrapolateElement(loc_asms[i], var, _local_to_global[i].rows,
                           _nodal_values, counts);
    }

    _nodal_values.componentwiseDivide(counts);
}

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
void
LocalLinearLeastSquaresExtrapolator<GlobalVector, VariableEnum, LocalAssembler>::
calculateResiduals(LocalAssemblers const& loc_asms, VariableEnum var)
{
    _residuals.resize(loc_asms.size());

    for (std::size_t i=0; i<loc_asms.size(); ++i)
    {
        _residuals[i] = calculateResiudalElement(
                    loc_asms[i], var, _local_to_global[i].rows,
                    _nodal_values
                    );
    }
}

} // namespace ProcessLib
