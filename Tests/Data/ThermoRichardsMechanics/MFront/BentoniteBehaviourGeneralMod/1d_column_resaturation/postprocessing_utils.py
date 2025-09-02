from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv


class IterPvd:
    def __init__(self, pvd):
        self._reader = pv.get_reader(pvd)
        self._time_step = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._time_step == len(self.time_values):
            raise StopIteration

        self._reader.set_active_time_point(self._time_step)
        mesh = self._reader.read()[0]
        self._time_step += 1
        return mesh

    @property
    def time_values(self):
        return self._reader.time_values


class IterVtus:
    def __init__(self, vtus, time_values):
        assert len(vtus) == len(time_values)
        self._vtus = vtus
        self.time_values = time_values
        self._time_step = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._time_step == len(self.time_values):
            raise StopIteration

        mesh = pv.read(self._vtus[self._time_step])
        self._time_step += 1
        return mesh


def extract_data(mesh_iterator, line_mesh=None):
    ts = np.asarray(mesh_iterator.time_values)
    profiles = []

    agg_names = ["min", "max", "mean", "std"]
    agg_fcts = [np.min, np.max, np.mean, np.std]
    aggregations = {agg_n: [] for agg_n in agg_names}

    for mesh in mesh_iterator:
        # line profiles
        if line_mesh is None:
            bounds = np.array(mesh.bounds).reshape((3, -1))
            center = np.mean(bounds, axis=1)
            line_z = [
                [center[0], center[1], bounds[2, 0]],
                [center[0], center[1], bounds[2, 1]],
            ]
            line_mesh = pv.Line(*line_z, resolution=1000)

        profile = line_mesh.sample(mesh)

        profiles.append(profile.point_data)

        # aggregations
        for agg_n, agg_f in zip(agg_names, agg_fcts, strict=True):
            agg_rec = {}
            for n in mesh.array_names:
                d = mesh.get_array(n)
                d_agg = agg_f(d, axis=0)
                if np.size(d_agg) == 1:
                    agg_rec[n] = d_agg
                else:
                    for comp, v in enumerate(d_agg):
                        agg_rec[f"{n}[{comp}]"] = v

            aggregations[agg_n].append(agg_rec)

    dfs_aggregation = {}
    for agg_n, agg_recs in aggregations.items():
        dfs_aggregation[agg_n] = pd.DataFrame.from_records(
            agg_recs, index=pd.Index(ts, name="time")
        )
        dfs_aggregation[agg_n]["time"] = ts

    return ts, line_mesh, dfs_aggregation, profiles


def plot_time_series(dfs_aggregation, quantities=None):
    plot_path = Path("out/plot_agg")
    plot_path.mkdir(exist_ok=True)

    if quantities is None:
        quantities = dfs_aggregation["min"].columns

    for qty_name in quantities:
        qty_series = dfs_aggregation["min"][qty_name]
        if qty_name == "time":
            continue
        n_comp = np.size(qty_series.iloc[0])
        assert n_comp == 1  # components separated during reading

        fig, ax = plt.subplots()

        for agg_n in ["min", "max", "mean"]:
            df_agg = dfs_aggregation[agg_n]
            df_agg.plot("time", qty_name, ax=ax, label=agg_n)

        ax.set_title(qty_name)
        fig.savefig(plot_path / f"time_series_{qty_name}.png")
        plt.close(fig)


def check_quantity(quantity, dfs_aggregation, checked_quantities, value, abstol):
    success = True

    for agg_n in ("min", "max"):
        df_agg = dfs_aggregation[agg_n]
        s = np.allclose(value, df_agg[quantity], atol=abstol, rtol=0)
        if not s:
            success = False
            print(f"Error: check failed: {quantity} != {value} +/i {abstol} ({agg_n})")

    if success:
        print(f"all({quantity} = {value} ± {abstol})")

    if success:
        checked_quantities.append(quantity)

    return success


def check_zero_and_constant_quantities(dfs_aggregation):
    cqs = []  # checked quantities
    dfs = dfs_aggregation
    s = True  # success

    s &= check_quantity("displacement[0]", dfs, cqs, 0, 1e-15)
    s &= check_quantity("displacement[1]", dfs, cqs, 0, 1e-15)

    s &= check_quantity("epsilon_ip[0]", dfs, cqs, 0, 1e-15)
    s &= check_quantity("epsilon_ip[1]", dfs, cqs, 0, 1e-15)
    s &= check_quantity("epsilon_ip[3]", dfs, cqs, 0, 1e-15)
    s &= check_quantity("epsilon_ip[4]", dfs, cqs, 0, 1e-15)
    s &= check_quantity("epsilon_ip[5]", dfs, cqs, 0, 1e-15)

    s &= check_quantity("intrinsic_permeability_ip[0]", dfs, cqs, 3e-20, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[1]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[2]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[3]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[4]", dfs, cqs, 3e-20, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[5]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[6]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[7]", dfs, cqs, 0, 1e-30)
    s &= check_quantity("intrinsic_permeability_ip[8]", dfs, cqs, 3e-20, 1e-30)

    s &= check_quantity("relative_permeability_ip", dfs, cqs, 1, 1e-15)

    s &= check_quantity("sigma_total_ip[3]", dfs, cqs, 0, 1e-7)
    s &= check_quantity("sigma_total_ip[4]", dfs, cqs, 0, 1e-7)
    s &= check_quantity("sigma_total_ip[5]", dfs, cqs, 0, 1e-7)

    s &= check_quantity("velocity_ip[0]", dfs, cqs, 0, 1e-19)
    s &= check_quantity("velocity_ip[1]", dfs, cqs, 0, 1e-19)

    s &= check_quantity("viscosity_ip", dfs, cqs, 0.0018, 1e-15)

    s &= check_quantity("temperature_interpolated", dfs, cqs, 298.15, 1e-13)

    s &= check_quantity("liquid_density_ip", dfs, cqs, 1000, 1e-10)

    s &= check_quantity("HeatFlowRate", dfs, cqs, 0, 1e-15)

    return s, cqs


# Cf. https://stackoverflow.com/a/26026189
def find_nearest_sorted(array, value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (
        idx == len(array) or np.abs(value - array[idx - 1]) < np.abs(value - array[idx])
    ):
        return idx - 1

    return idx


def select_profiles(ts, profiles, ts_selected):
    idcs = np.array([find_nearest_sorted(ts, t) for t in ts_selected])
    return ts[idcs], [profiles[i] for i in idcs]


def plot_profiles(ts, line_mesh, profiles, quantities, out_dir, profiles_ref=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    coords = line_mesh.points[:, 2]

    for qty in quantities:
        if qty.endswith("]"):
            assert qty[-3] == "["
            qty_name = qty[:-3]
            comp = int(qty[-2])
        else:
            qty_name = qty
            comp = None

        if profiles_ref is None:
            fig, ax = plt.subplots()
            axs = [ax]
        else:
            fig, axs = plt.subplots(2, sharex=True)
            ax = axs[0]

        for a in axs:
            a.set_prop_cycle(color=mpl.colormaps["rainbow"](np.linspace(0, 1, len(ts))))

        for ti, (t, prof) in enumerate(zip(ts, profiles, strict=True)):
            qty_values = prof[qty_name] if comp is None else prof[qty_name][:, comp]
            (h,) = ax.plot(
                coords,
                qty_values,
                label=f"$t$ = {t:.3g} s",
            )

            if profiles_ref is not None:
                prof_ref = profiles_ref[ti]
                qty_ref_values = (
                    prof_ref[qty_name] if comp is None else prof_ref[qty_name][:, comp]
                )
                ax.plot(
                    coords,
                    qty_ref_values,
                    color=h.get_color(),
                    label="ref" if ti == len(ts) - 1 else None,
                    ls="--",
                    linewidth=3,
                )

                axs[-1].plot(
                    coords, qty_values - qty_ref_values, label=f"$t$ = {t:.3g} s"
                )

        for a in axs:
            a.legend(loc="upper left", bbox_to_anchor=(1, 1))
        axs[-1].set_xlabel("$z$ / m")
        ax.set_title(qty)

        if profiles_ref is not None:
            axs[-1].set_ylabel("actual - reference")

        fig.set_size_inches(8, 5)
        fig.subplots_adjust(right=0.75)
        fig.savefig(out_dir / f"profile_{qty}.png", dpi=200)
        plt.close(fig)


def is_checked_via_another_quantity(all_quantities, q):
    for q2 in [
        q.replace("_avg", "_ip"),
        q + "_ip",
        q + "_interpolated",
        q.replace("[", "_ip["),
    ]:
        if q2 != q and q2 in all_quantities:
            break
    else:
        return False

    print(f"{q} is already checked via {q2}")
    return True
