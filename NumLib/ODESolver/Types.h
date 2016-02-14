#pragma once

namespace NumLib
{

//! \addtogroup ODESolver
//! @{

// TODO subject to change
using IndexType = std::size_t;

enum class NonlinearSolverTag : bool { Picard, Newton };

enum class ODESystemTag : char
{
    FirstOrderImplicitQuasilinear
};

//! @}

}
