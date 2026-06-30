import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys  # pyright: ignore[reportMissingImports]

out_dir = Path(os.environ.get("OGS_TESTRUNNER_OUT_DIR", "_out"))
out_dir.mkdir(parents=True, exist_ok=True)

last_timestep = 63072000  # Time at which the simulation ends

# 1. Define semi-empirical time-dependent logarithmic growth and decay formulation for the heat generation (Hanson et al., 2022)

# Constants to control the shape of the heat generation function
A = 4.88  # peak heat generation in W/m³
B = 50 * 24 * 60 * 60  # shape constant for the peak heat generation in seconds
C = 5000 * 24 * 60 * 60  # shape constant for the peak heat generation in seconds
D = 180 * 24 * 60 * 60  # decay rate factor in seconds


def time_dependent_source(t):
    return A * ((t / (B + t)) * (C / (C + t))) * np.exp(-np.sqrt(t / D))


# 2. Calculate heat generation as a function of expended-energy (Hanson et al. 2013)


def energy_expended(t):
    source_values = time_dependent_source(t)
    energy_vals = cumulative_trapezoid(source_values, t, initial=0)
    return source_values, energy_vals


heat_generation, energy_expended_values = energy_expended(
    np.linspace(0, 2 * 365 * 24 * 3600, 2 * 365 * 24)
)


# Snap function, to overcome precision issues of point coordinates
def snap_coords(sensor_coords, key_coords, tolerance=1e-3):
    xs, ys, _zs = sensor_coords
    xk, yk, _zk = key_coords

    return np.hypot(xs - xk, ys - yk) <= tolerance


# Integration points close to actual sensors
sensor_coords_list = [
    (10.369722119417698, 4.315470053837926, 0.0),
    (10.26536317669089, 6.484529946162075, 0.0),
]

sensor_df_list = [
    pd.DataFrame(
        {"time_s": [], "energy_expended_MJ/m3": [], "heat_generation_W/m3": []}
    )
    for i in range(len(sensor_coords_list))
]

energy_dict = {}
source_is_started = {}
last_value = {}
last_t = {}
timesteps = []


def calculate_source_term_for_layer(t, t_start_of_layer, time_of_construction, coords):
    if coords not in energy_dict:
        energy_dict[coords] = 0
        last_value[coords] = 0
        last_t[coords] = 0
        source_is_started[coords] = False

    if t < t_start_of_layer + time_of_construction:
        value = 0.0

    elif not source_is_started[coords]:
        value = time_dependent_source(
            t - t_start_of_layer - time_of_construction
        )  # Initialize source calculations with t-dependent source because initial last_value = 0
        last_value[coords] = value

    else:
        value = np.interp(
            energy_dict[coords], energy_expended_values, heat_generation, right=0
        )  # Use energy-expended formulation
        last_value[coords] = value

    # Output some stored data to a file after the calculation
    if (
        timesteps[-1] == last_timestep
    ):  # Overwrite the csv in every iteration in the last timestep
        for sensor_index, sensor_coords in enumerate(sensor_coords_list):
            if snap_coords(
                sensor_coords, coords
            ):  # only write if the integration point is the sensor point
                sensor_df_list[sensor_index].to_csv(
                    out_dir / f"energy_expended_sensor_index_{sensor_index}.csv",
                    index=False,
                )
    return value


def _init_timestep(dt, t):
    """Write the values at the beginning of a new timestep."""
    for key, hg_value in last_value.items():
        energy_dict[key] = dt * hg_value + energy_dict[key]
        if energy_dict[key] > 0.0:
            source_is_started[key] = True
    for sensor_index, sensor_coords in enumerate(sensor_coords_list):
        for key, hg_value in last_value.items():
            if snap_coords(sensor_coords, key):
                new_data_row = pd.Series(
                    {
                        "time_s": t,
                        "energy_expended_MJ/m3": energy_dict[key]
                        * 1e-6,  # convert from J to MJ
                        "heat_generation_W/m3": hg_value,
                    }
                )
                sensors = sensor_df_list[sensor_index]
                sensors.loc[len(sensors)] = new_data_row
                sensor_df_list[sensor_index] = sensors


# 3. Define OpenGeoSys Python Source Term for the landfill layer


class SourceTerm(OpenGeoSys.SourceTerm):
    def __init__(self):
        super().__init__()
        self._t_old_old = 0.0
        self._t_old = 0.0

    def _check_timestep(self, t):
        if self._t_old < t and t not in timesteps:
            timesteps.append(t)
            dt = self._t_old - self._t_old_old
            _init_timestep(dt, self._t_old)

            self._t_old_old = self._t_old
            self._t_old = t
            timesteps.append(self._t_old)

        elif t < self._t_old:
            self._t_old = t
            timesteps[-1] = t

    def getFlux(self, t, coords, _primary_vars):
        self._check_timestep(t)
        Jac = [0.0]

        x, y, z = coords
        # landfill base at y=4 and construction rate of 2e-7 m/s
        time_of_construction = (y - 4) / 2e-7
        t_start_of_layer = 0

        value = calculate_source_term_for_layer(
            t, t_start_of_layer, time_of_construction, (x, y, z)
        )
        return (value, Jac)


# 4. Instantiate source-term object referenced in the OpenGeoSys .prj file

source = SourceTerm()
