#include <cstdio>
#include <gtest/gtest.h>
#include "MathLib/Nonlinear/Root1D.h"
#include "MaterialsLib/Adsorption/DensityHauer.h"

using namespace Adsorption;


namespace
{

const double T = 303.15; // K
const double R = 8.314;  // J/mol/K
const double M = 0.018;  // kg/mol
const double phi = 0.4;
const double rho_SR0 = 1160.0; // kg/m^3
// const double k = 6e-3;   // s^-1

const double C0 = 0.0;
const double pV0 = 1e2;

std::unique_ptr<Reaction> ads(new Adsorption::DensityHauer);

double f(const double pV)
{
    const double C = ads->get_equilibrium_loading(pV, T, M);
    return (phi-1)/phi * (C - C0) * rho_SR0 * R*T/M + pV0 - pV;
}

double f2(const double pV)
{
    const double C = ads->get_equilibrium_loading(pV, T, M);
    // return (phi-1)/phi * (C - C0) * rho_SR0 * R*T/M + pV0 - pV;
    return (pV-pV0) * M/R/T * phi + (1-phi) * (C - C0) * rho_SR0;
}


} // anonymous namespace



template <typename T>
class RegulaFalsiTestZeolite : public ::testing::Test
{
};


using namespace MathLib::Nonlinear;
typedef ::testing::Types<Unmodified, Illinois, Pegasus, AndersonBjorck> RegulaFalsiTypes;
TYPED_TEST_CASE(RegulaFalsiTestZeolite, RegulaFalsiTypes);



TYPED_TEST(RegulaFalsiTestZeolite, Zeolite)
{
    auto rf = makeRegulaFalsi<TypeParam>(f, 1e-8, pV0);

    std::printf(" 0 -- x ~ %14.7g, range = %14.7g\n", rf.getResult(), rf.getRange());

    for (unsigned n=0; n<10; ++n)
    {
        rf.step(1);
        std::printf("%2i --  x ~ %14.7g, range = %14.7g\n", n+1, rf.getResult(), rf.getRange());
    }
}

TYPED_TEST(RegulaFalsiTestZeolite, Zeolite2)
{
    // auto rf = MakeRegulaFalsi<TypeParam>::create(f2, 1e-8, pV0);
    auto rf = makeRegulaFalsi<TypeParam>(f2, 1e-8, pV0);

    std::printf(" 0 -- x ~ %14.7g, range = %14.7g\n", rf.getResult(), rf.getRange());

    for (unsigned n=0; n<10; ++n)
    {
        rf.step(1);
        std::printf("%2i --  x ~ %14.7g, range = %14.7g\n", n+1, rf.getResult(), rf.getRange());
    }
}
