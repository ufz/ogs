#!/usr/bin/env python

import argparse
import logging
import os
import re
import sys
from collections import defaultdict
from decimal import ROUND_DOWN, ROUND_UP, Decimal
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

re_prj_file = re.compile("info: Reading project file (.*)[.]$")
re_line_compare = re.compile(
    "Comparing data array `(?P<array_1>.+)' from file `(?P<file_1>.+)' to data array `(?P<array_2>.+)' from file `(?P<file_2>.+)'[.]"
)
re_line_max_norm = re.compile(r"(abs|rel) maximum norm = \[([-+.,0-9e infa]+)\]")
re_line_error = re.compile(
    r"Absolute and relative error [(]maximum norm[)] are larger than the corresponding thresholds .* and .*[.]"
)

# A typical OGS output file name is, e.g., ramped_Neumann_BC_ts_50_t_0.625000.vtu.
# This regex matches the numerical parts of the file name including underscore
# before and underscore/dot after the number.
re_identify_number = re.compile(r"_(?P<num>[-+0-9.e]+)(?P<post>[_.])")

re_tests_data_dir = re.compile(r".*[\\/]Tests[\\/]Data[\\/]")


def _round_2_digits(d, rounding_mode):
    sign, dig, exp = d.as_tuple()
    d2 = Decimal(
        (0, dig, 1 - len(dig))
    )  # sign=0 is positive -> absolute value will be rounded
    d2_round = d2.quantize(Decimal("0.1"), rounding=rounding_mode)
    _sign2, dig2, exp2 = d2_round.as_tuple()
    d3 = Decimal((sign, dig2, exp2 - 1 + len(dig) + exp))
    # print(f"-+-> d={d}, {rounding_mode}\n +-> sign={sign}, dig={dig}, len(dig)={len(dig)}, exp={exp}\n +-> sign2={sign2}, dig2={dig2}, len(dig2)={len(dig2)}, exp2={exp2}\n +-> d3={float(d3)}")

    return float(d3)


def _round_up_2_digits(d):
    d_ = Decimal(d)

    if not d_.is_normal():
        return d

    # work around an error in Python's Decimal implementation
    up, _down = (ROUND_UP, ROUND_DOWN) if d >= 0 else (ROUND_DOWN, ROUND_UP)

    return _round_2_digits(d_, up)


def remove_trailing_whitespace_and_write_lines(fh, s):
    for line in s.splitlines():
        fh.write(line.rstrip() + "\n")


def parse_ogs_log(log_fh):
    try:
        log_fn = log_fh.name
    except AttributeError:
        log_fn = "<UNKNOWN>"

    mode = None
    vtkdiff_recs = []
    prj_file = None

    for line in log_fh:
        line = line.strip()  # noqa: PLW2901

        if mode is None:
            if line.endswith("info: ---------- vtkdiff begin ----------"):
                mode = "vtkdiff"
                vtkdiff_rec = {"check_succeeded": True}
            elif m := re_prj_file.search(line):
                prj_file = m.group(1).replace("\\", "/")
                logger.debug(
                    "The prj file that produced the log file '%s' is '%s'.",
                    log_fn,
                    prj_file,
                )

        elif mode == "vtkdiff":
            if line.endswith("info: ---------- vtkdiff end ----------"):
                mode = None
                vtkdiff_recs.append(vtkdiff_rec)
                del vtkdiff_rec
            elif m := re_line_compare.search(line):
                vtkdiff_rec["array_1"] = m.group("array_1")
                vtkdiff_rec["file_1"] = m.group("file_1").replace("\\", "/")
                vtkdiff_rec["array_2"] = m.group("array_2")
                vtkdiff_rec["file_2"] = m.group("file_2").replace("\\", "/")
            elif m := re_line_max_norm.search(line):
                abs_or_rel = m.group(1)
                error_norm = [float(comp) for comp in m.group(2).split(",")]
                vtkdiff_rec[abs_or_rel + "_maximum_norm"] = error_norm
            elif re_line_error.search(line):
                vtkdiff_rec["check_succeeded"] = False

    if mode is not None:
        logger.warning("File '%s' contains an incomplete vtkdiff block.", log_fn)

    logger.debug(
        "Log file '%s' contains %i vtkdiff records.", log_fn, len(vtkdiff_recs)
    )

    # split norm components
    recs = []
    for rec in vtkdiff_recs:
        try:
            abs_norm = rec["abs_maximum_norm"]
            rel_norm = rec["rel_maximum_norm"]
        except KeyError:
            continue
        assert len(abs_norm) == len(rel_norm)
        for comp, (a, r) in enumerate(zip(abs_norm, rel_norm, strict=True)):
            new_rec = dict(rec)
            new_rec["abs_maximum_norm"] = a
            new_rec["rel_maximum_norm"] = r
            new_rec["comp"] = comp
            recs.append(new_rec)

    return recs, prj_file


def agg_by_field(df):
    agg = {
        col: "min" if col == "check_succeeded" else "max"
        for col in df.columns
        if col != "array"
    }
    return df.groupby("array").agg(agg)


def round_up_2_digits(df):
    df2 = df.copy(deep=True)
    df2["abs_maximum_norm"] = df2["abs_maximum_norm"].map(_round_up_2_digits)
    df2["rel_maximum_norm"] = df2["rel_maximum_norm"].map(_round_up_2_digits)
    return df2


def remove_duplicate_columns(df):
    df2 = df.copy(deep=True)
    if (df2["array_1"] == df2["array_2"]).all():
        df2 = df2.drop(columns=["array_1", "array_2"])
        df2["array"] = df["array_1"]

    f1_base = df2["file_1"].map(os.path.basename)
    f2_base = df2["file_2"].map(os.path.basename)
    if (f1_base == f2_base).all():
        df2 = df2.drop(columns=["file_2"])
        df2 = df2.rename(columns={"file_1": "file"})
    return df2


def numbers_to_regexes(s):
    def replace(m):
        if m["post"] == "_":
            return "_.*_"
        return "_.*[.]"

    return re_identify_number.sub(replace, s)


def write_xml_snippet(df_max, xml_out_file):
    with xml_out_file.open("w") as fh:
        fh.write("<!-- summary of aggregated error norms:\n")
        print_table_summary(df_max, fh)

        if "file_1" in df_max.columns:
            assert "file_2" in df_max.columns
            fh.write("\n")
            fh.write(
                "    NOTE: File names of reference solution and actual solution differ!\n"
            )
            fh.write("          The simple regexes below might not work.\n")
            df_max = df_max.rename(columns={"file_1": "file"})

        fh.write("-->\n\n")

        fh.write("    <test_definition>")
        for tup in df_max.itertuples():
            file = Path(tup.file).name
            file_re = numbers_to_regexes(file)
            error_msg = (
                """
        <!-- Check failed -->"""
                if not tup.check_succeeded
                else ""
            )

            fh.write(
                f"""{error_msg}
        <vtkdiff>
            <regex>{file_re}</regex>
            <!-- <file>{file}</file> -->
            <field>{tup.Index}</field>
            <absolute_tolerance>{tup.abs_maximum_norm}</absolute_tolerance>
            <relative_tolerance>0.0 <!-- {tup.rel_maximum_norm} --></relative_tolerance>
        </vtkdiff>"""
            )
        fh.write("\n    </test_definition>\n")


def parse_ogs_log_files(path, *, decend_dirs=False):
    assert path.exists()

    if path.is_file():
        logger.debug("reading OGS log file '%s'.", path)
        vtkdiff_recs, prj_file = parse_ogs_log(path.open())
        if prj_file is None:
            return
        if not vtkdiff_recs:
            return
        yield (path.name, prj_file, pd.DataFrame.from_dict(vtkdiff_recs))
    elif decend_dirs and path.is_dir():
        logger.debug("looking for OGS log files in dir '%s'.", path)
        for path_ in path.iterdir():
            if path_.is_file():
                yield from parse_ogs_log_files(path_)


def group_by_prj_file(logs_prjs_dfs_vtkdiff):
    map_prj_file_to_logs_dfs_vtkdiff = defaultdict(list)

    for log_file, prj_file, df_vtkdiff in logs_prjs_dfs_vtkdiff:
        prj_rel_to_tests_data, num_subs = re_tests_data_dir.subn("", prj_file, count=1)
        if num_subs == 0:
            logger.warning(
                "The project file '%s' is not located under the Tests%sData directory.",
                prj_file,
                os.sep,
            )
            continue

        map_prj_file_to_logs_dfs_vtkdiff[prj_rel_to_tests_data].append(
            (log_file, df_vtkdiff)
        )

    return map_prj_file_to_logs_dfs_vtkdiff


def combine_all_vtkdiff_dataframes_for_each_prj_file(map_prj_file_to_logs_dfs_vtkdiff):
    res = {}

    for prj_file, logs_dfs_vtkdiff in map_prj_file_to_logs_dfs_vtkdiff.items():
        res[prj_file] = pd.concat((df for log, df in logs_dfs_vtkdiff))

    return res


def print_table_summary(df_max, fh):
    # we want neither file names nor components in the output
    df_for_output = df_max.drop(
        columns=["file", "file_1", "file_2", "comp"], errors="ignore"
    )

    with pd.option_context(
        "display.precision",
        1,
        "display.max_columns",
        10,
        "display.width",
        200,
        "display.max_colwidth",
        150,
    ) as _ctx:
        remove_trailing_whitespace_and_write_lines(fh, str(df_for_output))


def aggregate_log_files(log_files_dirs):
    # read log files
    logs_prjs_dfs_vtkdiff = [
        tup
        for log_file_dir in log_files_dirs
        for tup in parse_ogs_log_files(log_file_dir, decend_dirs=True)
    ]  # tuples (log file, prj file, aggregated vtkdiff records)

    map_prj_file_to_logs_dfs_vtkdiff = group_by_prj_file(logs_prjs_dfs_vtkdiff)

    map_prj_file_to_df_vtkdiff_combined = (
        combine_all_vtkdiff_dataframes_for_each_prj_file(
            map_prj_file_to_logs_dfs_vtkdiff
        )
    )

    map_prj_file_to_agg_vtkdiff_stats = {}

    for prj_file, df_vtkdiff in map_prj_file_to_df_vtkdiff_combined.items():
        df_vtkdiff = remove_duplicate_columns(df_vtkdiff)  # noqa: PLW2901

        map_prj_file_to_agg_vtkdiff_stats[prj_file] = agg_by_field(df_vtkdiff)

    return map_prj_file_to_agg_vtkdiff_stats


def run(log_files_dirs, out_dir, *, snippet_out, csv_out, verbose):
    map_prj_file_to_agg_vtkdiff_stats = aggregate_log_files(log_files_dirs)

    for i, (prj_file, df_vtkdiff_max) in enumerate(
        map_prj_file_to_agg_vtkdiff_stats.items()
    ):
        if verbose:
            if i != 0:
                print()  # blank line as a separator
            print(f"###### {prj_file}\n")

        df_vtkdiff_max_rounded = round_up_2_digits(df_vtkdiff_max)

        if verbose:
            print_table_summary(df_vtkdiff_max_rounded, sys.stdout)

        if out_dir is not None and (snippet_out or csv_out):
            prj_dir = Path(prj_file).parent
            assert not prj_dir.is_absolute()
            this_out_dir = out_dir / prj_dir
            this_out_dir.mkdir(parents=True, exist_ok=True)

            if snippet_out:
                this_xml_out_file = this_out_dir / Path(prj_file).name
                write_xml_snippet(df_vtkdiff_max_rounded, this_xml_out_file)

            if csv_out:
                this_csv_out_file = this_out_dir / f"{Path(prj_file).name}.csv"
                # for 17 digits precision cf. also https://en.cppreference.com/w/cpp/types/numeric_limits/max_digits10
                df_vtkdiff_max.to_csv(this_csv_out_file, sep="\t", float_format="%.17g")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Aggregates error norms in vtkdiff output in OGS log files and outputs XML snippets suitable as <test_definition>s in OGS prj files."
    )

    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Produce verbose cmdline output"
    )
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        help="output directory root for prj file snippets and CSV files.",
    )
    parser.add_argument(
        "--snippet-out",
        action="store_true",
        help="Output <test_definition> snippets for prj files.",
    )
    parser.add_argument(
        "--csv-out",
        action="store_true",
        help="Output aggregated statistics as CSV files.",
    )
    parser.add_argument(
        "log_files_dirs",
        type=Path,
        help="OGS log files or directories containing them",
        nargs="+",
    )

    args = parser.parse_args()

    snippet_out = args.snippet_out
    csv_out = args.csv_out
    outdir = args.out

    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if outdir and not (snippet_out or csv_out):
        logger.INFO("Output directory will not be used.")

    if (snippet_out or csv_out) and not outdir:
        logger.ERROR("No output directory specified for XML and/or CSV output.")
        sys.exit()

    if outdir is not None:
        assert outdir.is_dir()

    # suppress error message from Python interpreter, e.g., if a command
    # pipeline exits early, see
    # https://docs.python.org/3/library/signal.html#note-on-sigpipe
    try:
        run(
            args.log_files_dirs,
            outdir,
            snippet_out=snippet_out,
            csv_out=csv_out,
            verbose=args.verbose,
        )
    except BrokenPipeError:
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
