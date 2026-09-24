import numpy as np
import pandas as pd

# Total borehole/domain length (m); must equal LAYER_DZ.sum() (defined
# below, asserted there) -- the single source of truth for the layer
# geometry (see the LAYERS comment).
DOMAIN_LENGTH = 301


def add_MT_DB_geoboundaries_yaxis(
    ax,
    color="orange",
    style="--",
    alpha=0.25,
    lw=2,
    fs=10,
    tx=0,
    texts=True,
    masl=False,
):
    masl_factor = 1 / np.sqrt(2) if masl else 1
    # Hauptrogenstein
    ax.axhline(
        y=(DOMAIN_LENGTH - 37) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 264, 69
    # Passwang
    ax.axhline(
        y=(DOMAIN_LENGTH - 106) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 195, 32
    # Opa sandy I
    ax.axhline(
        y=(DOMAIN_LENGTH - 138) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 163, 36
    # Opa shaly I
    ax.axhline(
        y=(DOMAIN_LENGTH - 174) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 127, 12
    # Opa sandy II
    ax.axhline(
        y=(DOMAIN_LENGTH - 186) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 115, 4
    # Opa carbon
    ax.axhline(
        y=(DOMAIN_LENGTH - 190) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 111, 47
    # Opa shaly II
    ax.axhline(
        y=(DOMAIN_LENGTH - 237) * masl_factor,
        c=color,
        ls=style,
        alpha=alpha,
        lw=lw,
        zorder=99,
    )  # 64
    # Staffelegg
    # ax.axvline(x=-37, c=color, ls=style, alpha=alpha)
    if texts:
        ax.text(
            tx,
            (275) * masl_factor,
            "Hauptrogenstein",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (225) * masl_factor,
            "Passwang",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (175) * masl_factor,
            "Opalinus sandy I",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (145) * masl_factor,
            "Opalinus shaly I",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (120) * masl_factor,
            "Opalinus sandy II",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (110) * masl_factor,
            "Opalinus carbonate",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (85) * masl_factor,
            "Opalinus shaly II",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
        ax.text(
            tx,
            (25) * masl_factor,
            "Staffelegg",
            size=fs,
            rotation="horizontal",
            horizontalalignment="center",
        )
    return ax


def format_figure(axs):
    for idx in [0, 1]:
        add_MT_DB_geoboundaries_yaxis(
            axs[idx], color="red", style="--", tx=0.75e6, texts=False
        )
        axs_leg = axs[idx].legend(fontsize=10, frameon=True, ncols=1, loc="upper right")
        axs_leg.set_zorder(999999)
        axs_leg.get_frame().set_facecolor("white")
        axs_leg.get_frame().set_alpha(1.0)
        axs[idx].grid(visible=True, which="both", axis="both")

    axs[0].set_ylim(top=310)
    axs[0].set_title("A: Pressure")
    axs[0].set_ylabel("s / m")
    axs[0].set_xlabel("p / MPa")
    axs[1].set_title("B: Temperature")
    axs[1].set_xlabel("T / °C")
    axs[2].set_title("Model units")
    axs[2].set_xticks([])
    add_MT_DB_geoboundaries_yaxis(axs[2], color="red", style="--", tx=0.5, texts=True)


# Observation data
def load_observation_data_pressure():
    pressure_reference = pd.read_csv(
        "Observation_data/MTDB_modified_Goncalves_et_al_2023_pressure.csv"
    )
    p_ref_p = pressure_reference["X"].to_numpy() / 1e6
    p_ref_z = pressure_reference["Y_corr"].to_numpy()

    return p_ref_p, p_ref_z


def load_observation_data_temperature():
    temperature_fig_6 = pd.read_csv(
        "Observation_data/MTDB_modified_Goncalves_et_al_2023_temperature.csv"
    )
    t_ref_t_f6 = temperature_fig_6["X"].to_numpy()
    t_ref_z_f6 = temperature_fig_6["Y_corr"].to_numpy()
    return t_ref_t_f6, t_ref_z_f6


# Analytical solution for pressure
# Values from Gonçalvès 2023. All arrays from Staffelegg to Passwang, after flipping

gamma_w = 1e3 * 9.81  # N/m³
T_top = 12.43
ptop = 0.562e6
pbot = 2.463e6  # Pa
hbot = 530.0

dp = (pbot - ptop) - gamma_w * (DOMAIN_LENGTH / np.sqrt(2))

mu = 1e-3

# Layer thickness (dz, in m) and fitted temperature gradient (dT_dz, in
# K/m), Hauptrogenstein (top) to Staffelegg (bottom). Single source shared
# by analytical_solution_pressure's 'modified' variant, temperature_profile
# and _mesh.py's mesh geometry (_mesh.py imports LAYER_DZ), so a change to
# the layer geometry only has to be made in one place.
# See Fig. 6 in (Kiszkurno et al., 2026) for the fitted temperature
# gradients and Fig. 9 for gradient of intrinsic permeability.
LAYERS = [
    # name, dz, dT_dz
    ("Hauptrogenstein", 37, -0.05444633),
    ("Passwang I", 18, -0.05435378),
    ("Passwang II", 17, -0.05435378),
    ("Passwang III", 17, -0.05435378),
    ("Passwang IV", 17, -0.05435378),
    ("Opalinus sandy I", 32, -0.0489184),
    ("Opalinus shaly I", 36, -0.09381519),
    ("Opalinus sandy II", 12, -0.03280995),
    ("Opalinus carbonate", 4, -0.03280995),
    ("Opalinus shaly II", 47, -0.03280995),
    ("Staffelegg", 64, -0.06195801),
]
LAYER_DZ = np.array([dz for _, dz, _ in LAYERS])
LAYER_DT_DZ = np.array([dT_dz for _, _, dT_dz in LAYERS])
assert LAYER_DZ.sum() == DOMAIN_LENGTH


def _solve(kf, dz, eT, dT_dz, dh):
    n = len(kf)
    A = np.zeros((n + 1, n + 1))
    A[0:-1, 0:-1] += np.diag(kf / dz)
    A[0:-1, -1] += 1
    A[-1, 0:-1] += np.ones(n)
    RHS = np.zeros(n + 1)
    RHS[0:-1] = -eT * kf * dT_dz
    RHS[-1] = dh
    return np.linalg.solve(A, RHS)


def analytical_solution_pressure(base=True, TO=False, gradT=False):

    _, p_ref_z = load_observation_data_pressure()

    if gradT and base:
        dT_dz = np.array(
            [
                # See Fig. 6 in (Kiszkurno et al., 2026)
                # for the fitted temperature gradients
                -0.06195801,  # Staffelegg
                -0.03280995,  # Opalinus shaly II
                -0.03280995,  # Opalinus carbonate
                -0.03280995,  # Opalinus sandy II
                -0.09381519,  # Opalinus shaly I
                -0.0489184,  # Opalinus sandy I
                -0.05435378,  # Passwang
                -0.05435378,  # Hauptrogenstein
            ]
        )
    elif gradT and not base:
        # Staffelegg (top of the array) to Hauptrogenstein (bottom),
        # the reverse of LAYERS's Hauptrogenstein-to-Staffelegg order.
        dT_dz = np.flip(LAYER_DT_DZ)
    else:
        dT_dz = -0.0556
    if base:
        # Calculates base variant
        k = np.array(
            [
                # See Tab. 2 in (Kiszkurno et al.,  2026)
                # and Tab. 4 in (Gonçalvès et al., 2023)
                9.8e-18,  # Staffelegg
                7.8e-20,  # Opalinus shaly II
                9.06e-21,  # Opalinus carbonate
                3.74e-21,  # Opalinus sandy II
                4.76e-20,  # Opalinus shaly I
                3.74e-21,  # Opalinus sandy I
                1e-19,  # Passwang
                1e-10,  # Hauptrogenstein
            ]
        )  # m²
        eT = np.array(
            [
                # See Tab. 2 in (Kiszkurno et al.,  2026)
                # and Tab. 4 in (Gonçalvès et al., 2023)
                2.04e5,  # Staffelegg
                2.26e5,  # Opalinus shaly II
                7.78e4,  # Opalinus carbonate
                3.82e5,  # Opalinus sandy II
                1.53e5,  # Opalinus shaly I
                1.85e5,  # Opalinus sandy I
                0,  # Passwang
                0,  # Hauptrogenstein
            ]
        )  # Pa/K
        dz = np.array(
            [
                64,  # Staffelegg
                47,  # Opalinus shaly II
                4,  # Opalinus carbonate
                12,  # Opalinus sandy II
                36,  # Opalinus shaly I
                32,  # Opalinus sandy I
                69,  # Passwang
                37,  # Hauptrogenstein
            ]
        )
        kf = k / mu * gamma_w  # has to be in m/s
        eT /= gamma_w  # conversion to m/K
    else:
        # Calculates modified variant
        # With TO and permeability gradient in Passwang
        k = np.array(
            [
                # See Tab. 2 and 3 in (Kiszkurno et al.,  2026)
                # and Tab. 4 in (Gonçalvès et al., 2023)
                9.8e-18,  # Staffelegg
                7.8e-20,  # Opalinus shaly II
                9.06e-21,  # Opalinus carbonate
                3.74e-21,  # Opalinus sandy II
                4.76e-20,  # Opalinus shaly I
                3.74e-21,  # Opalinus sandy I
                1e-19,  # Passwang IV
                1e-20,  # Passwang III
                1e-21,  # Passwang II
                1e-20,  # Passwang I
                1e-10,  # Hauptrogenstein
            ]
        )
        eT = np.array(
            [
                # See Fig. 10 in (Kiszkurno et al., 2026)
                282865.75596207,  # Staffelegg
                177869.77885112,  # Opalinus shaly II
                38900.0,  # Opalinus carbonate
                191000.0,  # Opalinus sandy II
                76500.0,  # Opalinus shaly I
                169054.71157144,  # Opalinus sandy I
                277500.0,  # Passwang IV
                277500.0,  # Passwang III
                277500.0,  # Passwang II
                277500.0,  # Passwang I
                0.0,  # Hauptrogenstein
            ]
        )
        # Staffelegg (top of the array) to Hauptrogenstein (bottom),
        # the reverse of LAYERS's Hauptrogenstein-to-Staffelegg order.
        dz = np.flip(LAYER_DZ)
        kf = k / mu * gamma_w  # has to be in m/s
        eT /= gamma_w  # conversion to m/K
    if not TO:
        eT = np.zeros_like(eT)

    dh = -dp / gamma_w  # m
    solution = _solve(kf, dz, eT, dT_dz, dh)
    h_no_TO = hbot + np.append(0, solution[0:-1].cumsum())
    s = np.append(0, dz.cumsum())
    z = s / np.sqrt(2)
    p_sol = ((h_no_TO - hbot - z) * gamma_w + pbot) / 1e6
    z_sol = s

    z_out = DOMAIN_LENGTH + p_ref_z
    p_out = np.interp(z_out, z_sol, p_sol)

    return p_out, z_out


# Analytical solution for temperature


def temperature_profile(gradT=True):
    """Cumulative-gradient temperature profile shared by the mesh's initial
    condition (_mesh.py's _insert_T_init) and the benchmark's analytical
    curve (analytical_solution_temperature below); both reduce the same
    LAYERS dz/dT_dz data (module-level, Hauptrogenstein to Staffelegg) to a
    piecewise-linear T(z), only their consumers' post-processing differs.
    """

    if gradT:
        dT_dz = LAYER_DT_DZ
        dz = LAYER_DZ
    else:
        dz = np.array(
            [
                37,  # Hauptrogenstein
                69,  # Passwang
                32,  # Opalinus sandy I
                36,  # Opalinus shaly I
                12,  # Opalinus sandy II
                4,  # Opalinus carbonate
                47,  # Opalinus shaly II
                64,  # Staffelegg
            ]
        )
        dT_dz = -0.0556 * np.ones_like(dz)

    Tz = np.zeros((len(dz) + 1,))
    T_temp = T_top
    Tz[0] = T_temp
    for idx in range(len(dz)):
        T_temp = T_temp - dz[idx] * dT_dz[idx]
        Tz[idx + 1] = T_temp

    return Tz, DOMAIN_LENGTH - np.concatenate([[0], dz]).cumsum()


def analytical_solution_temperature(gradT=True):

    _, t_ref_z_f6 = load_observation_data_temperature()

    Tz, z_grid = temperature_profile(gradT)

    T_out = np.interp(
        np.flip(t_ref_z_f6),
        np.flip(z_grid),
        np.flip(Tz),
    )

    return np.flip(T_out), t_ref_z_f6


# Test against stored reference data


def test_against_reference_data(THM_p_mod, THM_T_mod, THM_TO_p_mod, THM_TO_T_mod):
    # A fixed absolute tolerance (as in np.testing.assert_almost_equal's
    # `decimal`) does not scale with magnitude: it is appropriate for the
    # O(10) degC temperatures but far too tight for the O(1e6) Pa pressures,
    # where it rejects sub-ULP-class solver noise (~1e-11 relative, from
    # linear-solver/BLAS differences across machines) as if it were a
    # regression. Use a relative tolerance instead, sized well above that
    # noise floor while still far tighter than any real physics change.
    rtol = 1e-6
    ref_data = np.load("Reference_data.npz")

    # With all modifications
    np.testing.assert_allclose(THM_p_mod, ref_data["THM_p_mod"], rtol=rtol)
    np.testing.assert_allclose(THM_T_mod, ref_data["THM_T_mod"], rtol=rtol)
    np.testing.assert_allclose(THM_TO_p_mod, ref_data["THM_TO_p_mod"], rtol=rtol)
    np.testing.assert_allclose(THM_TO_T_mod, ref_data["THM_TO_T_mod"], rtol=rtol)

    print("All tests were passed!")


def add_thermo_osmosis(project):

    kT = np.array(
        [
            # See Fig. 10 in (Kiszkurno et al., 2026)
            "0.0",  # Hauptrogenstein
            "2.78e-12",  # Passwang I
            "2.78e-13",  # Passwang II
            "2.78e-12",  # Passwang III
            "2.78e-11",  # Passwang IV
            "6.32e-13",  # Opalinus sandy I
            "3.64e-12",  # Opalinus shaly I
            "7.14e-13",  # Opalinus sandy II
            "3.52e-13",  # Opalinus carbonate
            "1.39e-11",  # Opalinus shaly II
            "2.77e-9",  # Staffelegg
        ]
    )

    for kT_id, kT_str in enumerate(kT):
        # kT_id+1 - indexing in prj file starts at 1
        project.add_block(
            blocktag="property",
            block_attrib=None,
            parent_xpath=f"./media/medium[{kT_id+1}]/properties",
            taglist=["name", "type", "value"],
            textlist=["thermal_osmosis_coefficient", "Constant", kT_str],
        )

    project.replace_text(
        value="MTDB_THM_TO_Modified",
        xpath="./time_loop/output/prefix",
    )

    return project
