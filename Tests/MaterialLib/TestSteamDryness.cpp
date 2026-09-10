// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <string>

#include "MaterialLib/MPL/Utils/SteamDryness.h"
#include "Tests/MaterialLib/LogCapture.h"

namespace
{
// Saturation enthalpies of water and steam at about 1 MPa, in J/kg.
constexpr double h_sat_liquid = 762.683e3;
constexpr double h_sat_vapour = 2777.12e3;

using Tests::LogCapture;
}  // namespace

TEST(MaterialLibSteamDryness, TwoPhaseIsTheEnthalpyFraction)
{
    using namespace MaterialPropertyLib;

    for (double dryness : {0., 0.25, 0.5, 0.75, 1.})
    {
        double const enthalpy =
            h_sat_liquid + dryness * (h_sat_vapour - h_sat_liquid);

        EXPECT_NEAR(dryness, steamDryness(enthalpy, h_sat_liquid, h_sat_vapour),
                    1e-15)
            << "dryness = " << dryness;
    }
}

// Subcooled liquid is a state the callers represent off the saturation line,
// so the cap only removes the negative vapour mass fraction.
TEST(MaterialLibSteamDryness, SubcooledLiquidHasNoVapour)
{
    using namespace MaterialPropertyLib;

    EXPECT_EQ(0., steamDryness(h_sat_liquid - 1., h_sat_liquid, h_sat_vapour));
    EXPECT_EQ(0., steamDryness(0., h_sat_liquid, h_sat_vapour));
    EXPECT_EQ(0., steamDryness(-1e6, h_sat_liquid, h_sat_vapour));
}

// Superheated steam is capped as well, because the drift-flux closure that
// consumes the dryness is only defined up to one.
TEST(MaterialLibSteamDryness, SuperheatedSteamIsFullyDry)
{
    using namespace MaterialPropertyLib;

    EXPECT_EQ(1., steamDryness(h_sat_vapour + 1., h_sat_liquid, h_sat_vapour));
    EXPECT_EQ(1., steamDryness(1e7, h_sat_liquid, h_sat_vapour));
}

// The superheat that the cap discards is not represented anywhere else, so the
// cap is contracted to report itself. A silent cap would leave a section
// evaluated at the saturation state without any trace of it in the log.
TEST(MaterialLibSteamDryness, SuperheatedSteamIsReported)
{
    using namespace MaterialPropertyLib;

    {
        LogCapture log;
        steamDryness(1e7, h_sat_liquid, h_sat_vapour);

        EXPECT_NE(std::string::npos, log.text().find("superheated"))
            << "log: " << log.text();
        EXPECT_NE(std::string::npos, log.text().find("1e+07"))
            << "log: " << log.text();
    }

    // Neither of the two states the mixture does represent is reported: the
    // two-phase range is the regular case, and the subcooled cap only removes
    // a negative vapour fraction from a state the callers re-evaluate.
    {
        LogCapture log;
        steamDryness(0.5 * (h_sat_liquid + h_sat_vapour), h_sat_liquid,
                     h_sat_vapour);
        steamDryness(h_sat_liquid - 1., h_sat_liquid, h_sat_vapour);

        EXPECT_EQ("", log.text());
    }
}

// At the critical point the two saturation enthalpies coincide and the ratio
// is not finite. It has to reach the closure, which rejects it, rather than
// being turned into a plausible looking dryness by either cap.
TEST(MaterialLibSteamDryness, NonFiniteRatioPassesThrough)
{
    using namespace MaterialPropertyLib;

    double const h_crit = 2084.26e3;

    EXPECT_TRUE(std::isnan(steamDryness(h_crit, h_crit, h_crit)));
    EXPECT_TRUE(std::isnan(steamDryness(
        std::numeric_limits<double>::quiet_NaN(), h_sat_liquid, h_sat_vapour)));
}
