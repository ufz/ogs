#ifndef PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
#define PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H

#include "Extrapolator.h"

namespace ProcessLib
{

template<typename GlobalVector, typename VariableEnum, typename LocalAssembler>
class LocalLinearLeastSquaresExtrapolator
        : public Extrapolator<GlobalVector, VariableEnum, LocalAssembler>
{
public:
    using LocalAssemblers = typename Extrapolator<GlobalVector, VariableEnum, LocalAssembler>::LocalAssemblers;

    explicit LocalLinearLeastSquaresExtrapolator(
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
    GlobalVector _residuals;
    AssemblerLib::LocalToGlobalIndexMap const& _local_to_global;
};

}

#include "LocalLinearLeastSquaresExtrapolator-impl.h"

#endif // PROCESSLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
