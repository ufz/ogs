#!/usr/bin/env python3

import os
import xml.etree.ElementTree as ET

import pandas as pd


def parse_and_filter_prj_file(file_path):
    path = []  # path in the XML element tree
    objs = []  # meta info associated with each level in the path
    for event, elem in ET.iterparse(file_path, ["start", "end"]):
        if event == "start":
            objs.append({})
            path.append(elem.tag)
        elif event == "end":
            p = "/".join(path)
            assert path.pop() == elem.tag
            obj = objs.pop()
            obj["path"] = p

            if p == "OpenGeoSysProject/processes/process/type":
                parent = objs[-1]
                parent["type"] = elem.text
            elif p == "OpenGeoSysProject/processes/process":
                yield obj  # this is interesting, yield it
            elif p.startswith("OpenGeoSysProject/media/"):
                if elem.tag in {"type", "name"}:
                    parent = objs[-1]
                    parent[elem.tag] = elem.text
                elif elem.tag in {
                    "phase",
                    # "component",
                    "property",
                }:
                    yield obj  # this is interesting, yield it


def parse_all_prj_files(datadir):
    records = []
    map_file_path_to_pcs_type = {}

    for _i, (root, _dirs, files) in enumerate(os.walk(datadir)):
        for f in files:
            if not f.endswith(".prj"):
                continue

            file_path = f"{root}/{f}"

            have_pcs_type = False
            have_media = False
            recs = []

            for meta in parse_and_filter_prj_file(file_path):
                if (
                    meta["path"] == "OpenGeoSysProject/processes/process"
                    and "type" in meta
                ):
                    have_pcs_type = True
                    map_file_path_to_pcs_type[file_path] = meta["type"]
                elif meta["path"].startswith("OpenGeoSysProject/media/"):
                    have_media = True
                    meta["file_path"] = file_path
                    recs.append(meta)

            # only "complete" prj files having both a process type and media
            # info defined are considered
            # TODO implement support for included prj file snippets
            if have_pcs_type and have_media:
                records += recs

    df_n_t_p = pd.DataFrame.from_records(records).set_index("file_path")
    df_pcst = pd.DataFrame.from_dict(
        map_file_path_to_pcs_type, orient="index", columns=["pcs_type"]
    )
    return (
        df_n_t_p.join(df_pcst)
        .drop_duplicates()
        .sort_values(["pcs_type", "path", "name"])
        .reset_index(drop=True)
    )

    # columns: Name, Type, ~xml xPath, ProCeSs Type


def main(datadir, docauxdir):
    df_n_t_p_pcst = parse_all_prj_files(datadir)

    df_n_t_p_pcst.to_json(
        os.path.join(docauxdir, "ctest-media-info.json"), orient="records", indent=2
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="""\
This script extracts information about used media properties from all ctest prj
files and puts them into the DocAux directory.\
"""
    )
    parser.add_argument("datadir", help="the Tests/Data directory OGS's source code")
    parser.add_argument("docauxdir", help="the DocAux directory OGS's build directory")
    args = parser.parse_args()

    main(args.datadir, args.docauxdir)
