#ifndef PROCESSLIB_GLOBAL_LLSQ_EXTRAPOLATOR_H
#define PROCESSLIB_GLOBAL_LLSQ_EXTRAPOLATOR_H

#include "Extrapolator.h"

namespace ProcessLib
{

template<typename GlobalMatrix, typename GlobalVector,
         typename VariableEnum, typename LocalAssembler>
class GlobalLinearLeastSquaresExtrapolator
        : public Extrapolator<GlobalVector, VariableEnum, LocalAssembler>
{
public:
    using LocalAssemblers = typename Extrapolator<GlobalVector, VariableEnum, LocalAssembler>::LocalAssemblers;

    explicit GlobalLinearLeastSquaresExtrapolator(
            AssemblerLib::LocalToGlobalIndexMap const& local_to_global)
        : _nodal_values(local_to_global.dofSize()),
          _local_to_global(local_to_global)
    {}


    void extrapolate(LocalAssemblers const& loc_asms, VariableEnum var) override;

    void calculateResiduals(LocalAssemblers const& loc_asms, VariableEnum var) override;

    GlobalVector const& getNodalValues() const override { return _nodal_values; }

    GlobalVector const& getElementResiduals() const override { return _residuals; }

private:
    GlobalVector _nodal_values;
    GlobalVector _integration_point_values;
    GlobalVector _residuals;
    AssemblerLib::LocalToGlobalIndexMap const& _local_to_global;
};

}

#include "GlobalLinearLeastSquaresExtrapolator-impl.h"

#endif // PROCESSLIB_GLOBAL_LLSQ_EXTRAPOLATOR_H
