import contextlib
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv

# ATM only for symmetric tensors. TODO extend
map_tensor_component_name_to_ogs_index = {
    "XX": 0,
    "YY": 1,
    "ZZ": 2,
    "XY": 3,
    "YZ": 4,
    "XZ": 5,
}

map_ogs_index_to_component_name = {
    3: {0: "X", 1: "Y", 2: "Z"},
    6: {index: name for name, index in map_tensor_component_name_to_ogs_index.items()},
    9: {index: "XYZ"[index // 3] + "XYZ"[index % 3] for index in range(9)},
}


map_mfront_name_to_ogs = {
    "Strain": "epsilon_ip",
    "Stress": "sigma_total_ip",
    "LiquidPressure": "pressure",
    "Saturation": "saturation_ip",
}

map_ogs_name_to_mfront = {ogs: mfront for mfront, ogs in map_mfront_name_to_ogs.items()}


def parse_mtest_res_header(f):
    multi_comp_re = re.compile(r"([0-9]+)th component of the (.*) [(](.*)[)]")

    names = []
    with Path(f).open() as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            col, name = line.strip().split(": ")

            if m := multi_comp_re.match(name):
                name = f"{m.group(3)}"
            names.append(name)
    return names


def read_mtest_results(mtest_output_file):
    names = parse_mtest_res_header(mtest_output_file)
    df_mtest = pd.read_csv(
        mtest_output_file, sep=" ", comment="#", header=None, names=names, dtype=float
    )
    df_mtest.index = df_mtest.index.rename("time_step")
    return df_mtest


def read_single_element_ogs_results(pvd_file):
    reader = pv.get_reader(pvd_file)

    recs = []
    ts = reader.time_values

    for ti, _t in enumerate(ts):
        reader.set_active_time_point(ti)
        mesh = reader.read()[0]
        assert mesh.n_cells == 1
        if ti == 0:
            N = mesh.n_points
        else:
            assert mesh.n_points == N

        rec = dict(mesh.point_data)
        assert set(rec.keys()).isdisjoint(mesh.cell_data.keys())
        rec.update(mesh.cell_data)
        assert set(rec.keys()).isdisjoint(mesh.field_data.keys())
        rec.update(mesh.field_data)

        recs.append(rec)

    assert len(recs) != 0

    shapes = {field: value.shape for field, value in recs[0].items()}

    # assert same number of components and tuples for all timesteps
    for rec in recs[1:]:
        for field, value in rec.items():
            assert value.shape == shapes[field]

    return ts, recs


def aggregate(rec, f):
    return {k: f(v, axis=0) for k, v in rec.items()}


def aggregate_ogs_results(recs):
    names = ("min", "max", "mean", "std")
    fs = (np.min, np.max, np.mean, np.std)

    return {
        n: [aggregate(rec, f) for rec in recs] for n, f in zip(names, fs, strict=True)
    }
    # recs_min = [ aggregate(rec, np.min) for rec in recs ]
    # print(recs_min)


def to_dataframes_with_mfront_names(recs_ogs_agg, ts_ogs):
    # suffix _mfront denotes MFront column naming scheme
    recs_ogs_agg_mfront = {
        agg: ogs_names_to_mfront(recs) for agg, recs in recs_ogs_agg.items()
    }
    dfs_ogs_agg_mfront = {
        agg: pd.DataFrame.from_records(recs, index=pd.Index(ts_ogs, name="time"))
        for agg, recs in recs_ogs_agg_mfront.items()
    }
    for df in dfs_ogs_agg_mfront.values():
        df["time"] = ts_ogs

    return dfs_ogs_agg_mfront


def ogs_names_to_mfront(recs_ogs):
    def gen_mfront_keys_values(rec_ogs):
        for field_ogs, v in rec_ogs.items():
            for comp_idx, comp_value in enumerate(np.atleast_1d(v)):
                # print(ogs_name_to_mfront(field_ogs, comp_idx, np.size(v)), comp_value)
                yield ogs_name_to_mfront(field_ogs, comp_idx, np.size(v)), comp_value

    recs_mfront = []
    for rec in recs_ogs:
        recs_mfront.append(dict(gen_mfront_keys_values(rec)))

    return recs_mfront


def report_and_plot_summary(df_res, file_prefix):
    df_agg = df_res.agg(["min", "max"])
    rec_max = df_agg.loc["max"]

    ser_range = (rec_max - df_agg.loc["min"]).drop("time")

    print(
        "\n## The following quantities did not change over the course of the simulation:"
    )
    for qty in ser_range[ser_range == 0].index.array:
        print(f"  {qty} (constant = {rec_max[qty]:g})")

    print("\n## Changing quantities are plotted...")
    changing_quantities = ser_range[ser_range != 0].index.array

    plot(df_res, quantities=changing_quantities, file_prefix=file_prefix)


def common_columns(*dfs):
    if not dfs:
        return set()

    return set(dfs[0].columns.to_numpy()).intersection(
        *(df.columns.to_numpy() for df in dfs[1:])
    )


def plot(
    *dfs, labels=None, quantities="union", file_prefix="", file_suffix="", diff=None
):
    assert len(dfs) >= 1

    if isinstance(quantities, str):
        if quantities == "intersection":
            quantities = common_columns(*dfs)
        elif quantities == "union":
            quantities = set(dfs[0].columns.to_numpy()).union(
                *(df.columns.to_numpy() for df in dfs[1:])
            )
        else:
            msg = f"Unknown magic value for quantities argument: '{quantities}'"
            raise ValueError(msg)

    if diff is None:
        diff = len(dfs) > 1

    for qty in sorted(quantities):
        if qty == "time":
            continue

        if diff:
            fig, (ax, ax_diff) = plt.subplots(2, sharex=True)
            ax_last = ax_diff
            ax_diff.set_ylabel("absolute difference")

        else:
            fig, ax = plt.subplots()
            ax_last = ax
        # ax.set_ylabel(qty)

        for i, df in enumerate(dfs):
            assert "time" in df.columns
            try:
                kwargs = {}
                try:
                    kwargs["label"] = labels[i]
                except TypeError:
                    if len(dfs) > 1:
                        kwargs["label"] = f"df #{i}"

                df.plot("time", qty, ax=ax, **kwargs)
            except KeyError as e:
                print(f"KeyError in dfs[{i}]:", e)

        if diff:
            ax_diff.plot([], [])  # plot nothing to skip first plot style
            df_ref = dfs[0]
            ax_abs_max_095 = 0
            ax_abs_max = 0
            some_below_0 = False
            some_above_0 = False

            for i, df in enumerate(dfs[1:]):
                ts = df["time"].to_numpy()
                try:
                    assert np.allclose(
                        ts, df_ref["time"].to_numpy(), atol=1e-13, rtol=0
                    )
                except AssertionError:
                    print(ts - df_ref["time"].to_numpy())
                    raise

                kwargs = {}
                with contextlib.suppress(TypeError):
                    kwargs["label"] = f"{labels[i+1]} $-$ {labels[0]}"

                ds = df[qty].to_numpy() - df_ref[qty].to_numpy()

                ax_diff.plot(ts, ds, **kwargs)

                ds_abs = np.abs(ds)
                ax_abs_max_095 = max(ax_abs_max_095, np.quantile(ds_abs, 0.95))
                ax_abs_max = max(ax_abs_max, np.max(ds_abs))
                some_below_0 = some_below_0 or np.any(ds < 0)
                some_above_0 = some_above_0 or np.any(ds > 0)

            ax_diff.legend()
            if ax_abs_max != 0:
                if ax_abs_max_095 != 0:
                    linthresh = ax_abs_max_095
                else:
                    linthresh = 1e-15 * ax_abs_max

                ax_diff.set_yscale("symlog", linthresh=linthresh)
                if some_above_0:
                    ax_diff.axhline(linthresh, ls=":", color="gray")
                if some_below_0:
                    ax_diff.axhline(-linthresh, ls=":", color="gray")
            else:
                ax_diff.annotate(
                    "all differences are zero", (0.05, 0.05), xycoords="axes fraction"
                )
            ax_diff.axhline(0, ls=":", color="gray")

        ax.set_title(qty)
        ax_last.set_xlabel("$t$ / s")
        fig = ax.get_figure()
        print(f"  {file_prefix}{qty}{file_suffix}.png")
        fig.savefig(f"{file_prefix}{qty}{file_suffix}.png")
        plt.close(fig)


def mfront_to_ogs(rec):
    ogs_rec = {}

    for k, v in rec.items():
        if k == "time_step":
            continue

        name, comp = mfront_name_to_ogs(k)

        if comp is None:
            ogs_rec[name] = np.array([v])
        elif name not in ogs_rec:
            ogs_rec[name] = np.array([0] * (comp - 1) + [v])
        else:
            arr = ogs_rec[name]
            if len(arr) <= comp:
                arr.resize(comp + 1, refcheck=False)
            arr[comp] = v

    return ogs_rec


def mfront_name_to_ogs(k):
    assert k != "time"
    if k[-2:] in map_tensor_component_name_to_ogs_index:
        name = k[:-2]
        comp = map_tensor_component_name_to_ogs_index[k[-2:]]
    else:
        name = k
        comp = None

    if name in map_mfront_name_to_ogs:
        name = map_mfront_name_to_ogs[name]
    else:
        name = f"material_state_variable_{name}_ip"

    return name, comp


def ogs_name_to_mfront(name, comp, num_comps):
    if name.startswith("material_state_variable_") and name.endswith("_ip"):
        assert num_comps == 1  # other cases not yet implemented
        return name[len("material_state_variable_") : -len("_ip")]

    if num_comps == 1:
        comp_suffix = ""
    else:
        comp_suffix = map_ogs_index_to_component_name[num_comps][comp]

    name_mfront = map_ogs_name_to_mfront.get(name, name)

    return name_mfront + comp_suffix


def generate_reference_meshes(
    input_mesh, df_res_mtest, time_steps, filename_template, num_int_pts
):
    mesh_ic = pv.read(input_mesh)
    mesh_ic.point_data.clear()
    mesh_ic.cell_data.clear()
    mesh_ic.field_data.clear()
    N = mesh_ic.n_points

    df_res_t = df_res_mtest.copy(deep=False)
    df_res_t["time_step"] = df_res_t.index
    df_res_t = df_res_t.set_index("time")

    for t in time_steps:
        mesh = mesh_ic.copy(deep=False)
        rec = df_res_t.loc[t]
        for k, v in mfront_to_ogs(rec).items():
            if k.endswith("_ip"):
                mesh.field_data[k] = np.tile(v, (num_int_pts, 1))
            else:
                mesh.point_data[k] = np.tile(v, (N, 1))

        mesh.save(filename_template.format(t=t, ts=int(rec["time_step"])))


def write_test_definition_snippet(
    fields_mfront, fields_ogs_mfront, regex_prefix, output_file
):
    fields1 = {
        mfront_name_to_ogs(field)[0] for field in fields_mfront if field != "time"
    }
    fields2 = {
        mfront_name_to_ogs(field)[0] for field in fields_ogs_mfront if field != "time"
    }

    fields_ogs = fields1 & fields2

    with Path(output_file).open("w") as fh:
        fh.write("    <test_definition>")

        for field in sorted(fields_ogs):
            fh.write(
                f"""
        <vtkdiff>
            <regex>{regex_prefix}_ts_.*_t_.*[.]vtu</regex>
            <field>{field}</field>
            <absolute_tolerance>0.0</absolute_tolerance>
            <relative_tolerance>0.0</relative_tolerance>
        </vtkdiff>"""
            )

        fh.write("\n    </test_definition>\n")
