#include "LinearElasticIsotropicPhaseField.h"

namespace MaterialLib::Solids::Phasefield
{

/**
 * Double–minor symmetrized tensor product
 *
 * aOdotB computes the symmetric 4th-order tensor product
 * \f[
 * (A \odot B)_{ijkl}
 *   = \tfrac{1}{4}
 *     (A_{ik}B_{jl} + A_{il}B_{jk} + A_{jk}B_{il} + A_{jl}B_{ik}),
 * \f]
 * which is symmetric in both index pairs
 * \f$(i,j)\f$ and \f$(k,l)\f$ (“double-minor symmetry”).
 *
 * Reference
 * Mehrabadi, M.M., & Cowin, S.C. (1990). *Eigentensors of linear anisotropic
 * elastic materials.* Q. J. Mech. Appl. Math., 43(1), 15–41.
 */

template <>
MathLib::KelvinVector::KelvinMatrixType<3> aOdotB<3>(
    MathLib::KelvinVector::KelvinVectorType<3> const& A,
    MathLib::KelvinVector::KelvinVectorType<3> const& B)
{
    MathLib::KelvinVector::KelvinMatrixType<3> result;

    result(0, 0) = A(0) * B(0);
    result(0, 1) = result(1, 0) = A(3) * B(3) / 2.;
    result(0, 2) = result(2, 0) = A(5) * B(5) / 2.;
    result(0, 3) = result(3, 0) = (A(0) * B(3) + A(3) * B(0)) / 2;
    result(0, 4) = result(4, 0) =
        (A(3) * B(5) + A(5) * B(3)) / (2. * std::sqrt(2.));
    result(0, 5) = result(5, 0) = (A(0) * B(5) + A(5) * B(0)) / 2;

    result(1, 1) = A(1) * B(1);
    result(1, 2) = result(2, 1) = A(4) * B(4) / 2.;
    result(1, 3) = result(3, 1) = (A(1) * B(3) + A(3) * B(1)) / 2.;
    result(1, 4) = result(4, 1) = (A(1) * B(4) + A(4) * B(1)) / 2.;
    result(1, 5) = result(5, 1) =
        (A(3) * B(4) + A(4) * B(3)) / (2. * std::sqrt(2.));

    result(2, 2) = A(2) * B(2);
    result(2, 3) = result(3, 2) =
        (A(4) * B(5) + A(5) * B(4)) / (2. * std::sqrt(2.));
    result(2, 4) = result(4, 2) = (A(2) * B(4) + A(4) * B(2)) / 2.;
    result(2, 5) = result(5, 2) = (A(2) * B(5) + A(5) * B(2)) / 2.;

    result(3, 3) = (A(0) * B(1) + A(1) * B(0) + A(3) * B(3)) / 2.;
    result(3, 4) = result(4, 3) =
        (A(3) * B(4) + A(4) * B(3)) / 4. +
        (A(1) * B(5) + A(5) * B(1)) / (2. * std::sqrt(2.));
    result(3, 5) = result(5, 3) =
        (A(3) * B(5) + A(5) * B(3)) / 4. +
        (A(0) * B(4) + A(4) * B(0)) / (2. * std::sqrt(2.));

    result(4, 4) = (A(1) * B(2) + A(2) * B(1) + A(4) * B(4)) / 2.;
    result(4, 5) = result(5, 4) =
        (A(4) * B(5) + A(5) * B(4)) / 4. +
        (A(2) * B(3) + A(3) * B(2)) / (2. * std::sqrt(2.));

    result(5, 5) = (A(0) * B(2) + A(2) * B(0) + A(5) * B(5)) / 2.;
    return result;
}

template <>
MathLib::KelvinVector::KelvinMatrixType<2> aOdotB<2>(
    MathLib::KelvinVector::KelvinVectorType<2> const& A,
    MathLib::KelvinVector::KelvinVectorType<2> const& B)
{
    MathLib::KelvinVector::KelvinMatrixType<2> result;

    result(0, 0) = A(0) * B(0);
    result(0, 1) = result(1, 0) = A(3) * B(3) / 2.;
    result(0, 2) = result(2, 0) = 0;
    result(0, 3) = result(3, 0) = (A(0) * B(3) + A(3) * B(0)) / 2;

    result(1, 1) = A(1) * B(1);
    result(1, 2) = result(2, 1) = 0;
    result(1, 3) = result(3, 1) = (A(1) * B(3) + A(3) * B(1)) / 2.;

    result(2, 2) = A(2) * B(2);
    result(2, 3) = result(3, 2) = 0;

    result(3, 3) = (A(0) * B(1) + A(1) * B(0) + A(3) * B(3)) / 2.;

    return result;
}

}  // namespace MaterialLib::Solids::Phasefield
