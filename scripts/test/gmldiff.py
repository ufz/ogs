#!/usr/bin/env python3

# Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
#            Distributed under a Modified BSD License.
#              See accompanying file LICENSE.txt or
#              http://www.opengeosys.org/project/license

from xml.dom import minidom
import argparse
import math

parser = argparse.ArgumentParser(description="Diff OpenGeoSys GML files.")
parser.add_argument(
    "gmls", metavar="gml", type=str, nargs=2, help="Two gml files to compare"
)
parser.add_argument("--abs", type=float, default=1e-16)
parser.add_argument("--rel", type=float, default=1e-16)

args = parser.parse_args()

docA = minidom.parse(args.gmls[0])
docB = minidom.parse(args.gmls[1])

name = docA.getElementsByTagName("name")[0]
print(
    f"Comparing gml with name '{name.firstChild.data}', abs={args.abs}, rel={args.rel}"
)

pointsA = docA.getElementsByTagName("point")
pointsB = docB.getElementsByTagName("point")

if len(pointsA) != len(pointsB):
    print("Mismatch of number of points!")
    exit(1)

for pointA, pointB in zip(pointsA, pointsB):
    if int(pointA.getAttribute("id")) != int(pointB.getAttribute("id")):
        print("Points do not have the same order!")
        exit(1)

    print(pointB.getAttribute("id"))
    for dim in ["x", "y", "z"]:
        a = float(pointA.getAttribute(dim))
        b = float(pointB.getAttribute(dim))
        if not math.isclose(
            a,
            b,
            rel_tol=args.rel,
            abs_tol=args.abs,
        ):
            print(
                f"Point with id={pointA.getAttribute('id')} differ: abs={abs(a - b)}, rel={abs(a - b) / b}"
            )
            exit(1)

exit(0)
