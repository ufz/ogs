import math
from pathlib import Path

import gmsh


def mesh_GreatCell_intact(
    lc=0.0025,
    lc2=0.0025,
    r0=0.097,
    r1=0.094,
    r2=0.09,
    r3=0.065,
    fracture_angle_deg=0.0,
    out_dir=".",
    meshname="Greatcell_mesh",
    mode="domain",
):
    out_path = Path(out_dir)
    if not out_path.exists():
        out_path.mkdir(parents=True)
        print(f"Created output directory: {out_dir}")
    else:
        print(f"Output directory already exists: {out_dir}")

    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.model.add("circular_fracture")

    th = math.radians(fracture_angle_deg)
    sth = math.sin(th)
    cth = math.cos(th)

    def p(x, y, mesh_size, tag):
        gmsh.model.geo.addPoint(x, y, 0, mesh_size, tag)

    def circle(p1, pc, p2, tag):
        gmsh.model.geo.addCircleArc(p1, pc, p2, tag)

    # Outer points
    outer_coords = [
        (1, 0.0989400395497483, -0.0145213144685405),
        (2, 0.0989400395497483, 0.0),
        (3, 0.0989400395497483, 0.0145213144685405),
        (4, 0.0979528917494459, 0.0194840415895657),
        (5, 0.0969657439491436, 0.0244467687105909),
        (6, 0.0914086774858697, 0.0378627139332354),
        (7, 0.0858516110225958, 0.05127865915588),
        (8, 0.0830404519257676, 0.0554858560599625),
        (9, 0.0802292928289395, 0.0596930529640449),
        (10, 0.0699611728964922, 0.0699611728964922),
        (11, 0.0596930529640449, 0.0802292928289395),
        (12, 0.0554858560599625, 0.0830404519257676),
        (13, 0.05127865915588, 0.0858516110225958),
        (14, 0.0378627139332354, 0.0914086774858697),
        (15, 0.0244467687105908, 0.0969657439491436),
        (16, 0.0194840415895657, 0.0979528917494459),
        (17, 0.0145213144685405, 0.0989400395497483),
        (18, 0.0, 0.0989400395497483),
        (19, -0.0145213144685405, 0.0989400395497483),
        (20, -0.0194840415895657, 0.0979528917494459),
        (21, -0.0244467687105908, 0.0969657439491436),
        (22, -0.0378627139332354, 0.0914086774858697),
        (23, -0.05127865915588, 0.0858516110225958),
        (24, -0.0554858560599625, 0.0830404519257676),
        (25, -0.0596930529640449, 0.0802292928289395),
        (26, -0.0699611728964922, 0.0699611728964922),
        (27, -0.0802292928289395, 0.0596930529640449),
        (28, -0.0830404519257677, 0.0554858560599624),
        (29, -0.0858516110225958, 0.05127865915588),
        (30, -0.0914086774858697, 0.0378627139332354),
        (31, -0.0969657439491436, 0.0244467687105908),
        (32, -0.0979528917494459, 0.0194840415895657),
        (33, -0.0989400395497483, 0.0145213144685405),
        (34, -0.0989400395497483, 0.0),
        (35, -0.0989400395497483, -0.0145213144685405),
        (36, -0.0979528917494459, -0.0194840415895657),
        (37, -0.0969657439491436, -0.0244467687105909),
        (38, -0.0914086774858697, -0.0378627139332354),
        (39, -0.0858516110225958, -0.05127865915588),
        (40, -0.0830404519257676, -0.0554858560599625),
        (41, -0.0802292928289395, -0.0596930529640449),
        (42, -0.0699611728964922, -0.0699611728964922),
        (43, -0.059693052964045, -0.0802292928289395),
        (44, -0.0554858560599625, -0.0830404519257676),
        (45, -0.05127865915588, -0.0858516110225958),
        (46, -0.0378627139332354, -0.0914086774858697),
        (47, -0.0244467687105908, -0.0969657439491436),
        (48, -0.0194840415895657, -0.0979528917494459),
        (49, -0.0145213144685405, -0.0989400395497483),
        (50, 0.0, -0.0989400395497483),
        (51, 0.0145213144685405, -0.0989400395497483),
        (52, 0.0194840415895657, -0.0979528917494459),
        (53, 0.0244467687105909, -0.0969657439491436),
        (54, 0.0378627139332354, -0.0914086774858697),
        (55, 0.05127865915588, -0.0858516110225958),
        (56, 0.0554858560599625, -0.0830404519257676),
        (57, 0.0596930529640449, -0.0802292928289395),
        (58, 0.0699611728964922, -0.0699611728964922),
        (59, 0.0802292928289396, -0.0596930529640449),
        (60, 0.0830404519257677, -0.0554858560599624),
        (61, 0.0858516110225958, -0.0512786591558799),
        (62, 0.0914086774858697, -0.0378627139332354),
        (63, 0.0969657439491436, -0.0244467687105908),
        (64, 0.0979528917494459, -0.0194840415895657),
    ]
    for tag, x, y in outer_coords:
        p(x, y, lc, tag)
    for i in range(1, 65):
        gmsh.model.geo.addLine(i, i + 1 if i < 64 else 1, 100 + i)

    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 1000)

    # Center
    p(0, 0, lc2, 2000)

    rid = {"r0": r0, "r1": r1, "r2": r2, "r3": r3}
    circle_loops = {}
    offset = 3000

    for i, (rk, rad) in enumerate(sorted(rid.items(), key=lambda item: item[1])):
        start = offset + i * 10
        tags = {"E": start + 1, "N": start + 2, "W": start + 3, "S": start + 4}

        p(rad * cth, rad * sth, lc2, tags["E"])
        p(-rad * sth, rad * cth, lc2, tags["N"])
        p(-rad * cth, -rad * sth, lc2, tags["W"])
        p(rad * sth, -rad * cth, lc2, tags["S"])

        arcs = [start + 5, start + 6, start + 7, start + 8]
        circle(tags["E"], 2000, tags["N"], arcs[0])
        circle(tags["N"], 2000, tags["W"], arcs[1])
        circle(tags["W"], 2000, tags["S"], arcs[2])
        circle(tags["S"], 2000, tags["E"], arcs[3])

        loop_id = 1100 + i
        gmsh.model.geo.addCurveLoop(arcs, loop_id)
        circle_loops[rk] = loop_id

    boundary_lines = {
        "PEE1": [101, 102],
        "PEE2": [162, 161],
        "PEE3": [158, 157],
        "PEE4": [154, 153],
        "PEE5": [150, 149],
        "PEE6": [146, 145],
        "PEE7": [142, 141],
        "PEE8": [138, 137],
        "PEE1a": [134, 133],
        "PEE2a": [130, 129],
        "PEE3a": [126, 125],
        "PEE4a": [122, 121],
        "PEE5a": [118, 117],
        "PEE6a": [114, 113],
        "PEE7a": [110, 109],
        "PEE8a": [106, 105],
        "DSS1": [163, 164],
        "DSS2": [160, 159],
        "DSS3": [156, 155],
        "DSS4": [152, 151],
        "DSS5": [148, 147],
        "DSS6": [144, 143],
        "DSS7": [140, 139],
        "DSS8": [136, 135],
        "DSS1a": [132, 131],
        "DSS2a": [128, 127],
        "DSS3a": [124, 123],
        "DSS4a": [120, 119],
        "DSS5a": [116, 115],
        "DSS6a": [112, 111],
        "DSS7a": [108, 107],
        "DSS8a": [104, 103],
    }

    surface_id = 200
    radii_sorted = sorted(rid, key=rid.get)
    gmsh.model.geo.addPlaneSurface([circle_loops[radii_sorted[0]]], surface_id)

    for i in range(1, len(radii_sorted)):
        surface_id += 1
        gmsh.model.geo.addPlaneSurface(
            [circle_loops[radii_sorted[i]], circle_loops[radii_sorted[i - 1]]],
            surface_id,
        )

    # Outer boundary
    surface_id += 1
    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 3000)
    gmsh.model.geo.addPlaneSurface([3000, circle_loops[radii_sorted[-1]]], surface_id)

    # Recombine surfaces
    for sid in range(200, 204):
        gmsh.model.geo.mesh.setRecombine(2, sid)

    p(0.04, 0.0, lc2, 700)
    gmsh.model.geo.addLine(2000, 700, 1234)  # Fracture line

    gmsh.model.geo.synchronize()

    # Physical Embed fracture line
    gmsh.model.mesh.embed(1, [1234], 2, 200)
    # Physical Source point at center
    gmsh.model.mesh.embed(0, [2000], 2, 200)
    gmsh.model.mesh.embed(0, [700], 2, 200)
    gmsh.model.mesh.generate(2)

    if mode == "BC":
        for name, lines in boundary_lines.items():
            gmsh.model.addPhysicalGroup(1, lines, name=name)
        gmsh.model.addPhysicalGroup(0, [2000], name="Source")
        gmsh.model.addPhysicalGroup(0, [2], name="p_right")
        gmsh.model.addPhysicalGroup(0, [18], name="p_top")
        gmsh.model.addPhysicalGroup(0, [34], name="p_left")
        gmsh.model.addPhysicalGroup(0, [50], name="p_bottom")
    elif mode == "domain":
        # Physical surfaces
        gmsh.model.addPhysicalGroup(2, [200], name="Central_sample")
        gmsh.model.addPhysicalGroup(2, [201, 202, 203], name="OuterPart_sample")
        gmsh.model.addPhysicalGroup(2, [204], name="Rubber_sheath")

    msh_file = Path(out_path, f"{meshname}.msh")
    gmsh.write(str(msh_file))
    gmsh.finalize()
    return msh_file


def mesh_GreatCell_embeddedFracture(
    lc=0.0025,
    lc2=0.0025,
    r0=0.097,
    r1=0.094,
    r2=0.09,
    r3=0.065,
    fracture_angle_deg=0.0,
    out_dir=".",
    meshname="Greatcell_mesh",
    mode="domain",
):
    out_path = Path(out_dir)
    if not out_path.exists():
        out_path.mkdir(parents=True)
        print(f"Created output directory: {out_dir}")
    else:
        print(f"Output directory already exists: {out_dir}")

    if not gmsh.isInitialized():
        gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.model.add("circular_fracture")

    th = math.radians(fracture_angle_deg)
    sth = math.sin(th)
    cth = math.cos(th)

    def p(x, y, mesh_size, tag):
        gmsh.model.geo.addPoint(x, y, 0, mesh_size, tag)

    def circle(p1, pc, p2, tag):
        gmsh.model.geo.addCircleArc(p1, pc, p2, tag)

    # Outer points
    outer_coords = [
        (1, 0.0989400395497483, -0.0145213144685405),
        (2, 0.0989400395497483, 0.0),
        (3, 0.0989400395497483, 0.0145213144685405),
        (4, 0.0979528917494459, 0.0194840415895657),
        (5, 0.0969657439491436, 0.0244467687105909),
        (6, 0.0914086774858697, 0.0378627139332354),
        (7, 0.0858516110225958, 0.05127865915588),
        (8, 0.0830404519257676, 0.0554858560599625),
        (9, 0.0802292928289395, 0.0596930529640449),
        (10, 0.0699611728964922, 0.0699611728964922),
        (11, 0.0596930529640449, 0.0802292928289395),
        (12, 0.0554858560599625, 0.0830404519257676),
        (13, 0.05127865915588, 0.0858516110225958),
        (14, 0.0378627139332354, 0.0914086774858697),
        (15, 0.0244467687105908, 0.0969657439491436),
        (16, 0.0194840415895657, 0.0979528917494459),
        (17, 0.0145213144685405, 0.0989400395497483),
        (18, 0.0, 0.0989400395497483),
        (19, -0.0145213144685405, 0.0989400395497483),
        (20, -0.0194840415895657, 0.0979528917494459),
        (21, -0.0244467687105908, 0.0969657439491436),
        (22, -0.0378627139332354, 0.0914086774858697),
        (23, -0.05127865915588, 0.0858516110225958),
        (24, -0.0554858560599625, 0.0830404519257676),
        (25, -0.0596930529640449, 0.0802292928289395),
        (26, -0.0699611728964922, 0.0699611728964922),
        (27, -0.0802292928289395, 0.0596930529640449),
        (28, -0.0830404519257677, 0.0554858560599624),
        (29, -0.0858516110225958, 0.05127865915588),
        (30, -0.0914086774858697, 0.0378627139332354),
        (31, -0.0969657439491436, 0.0244467687105908),
        (32, -0.0979528917494459, 0.0194840415895657),
        (33, -0.0989400395497483, 0.0145213144685405),
        (34, -0.0989400395497483, 0.0),
        (35, -0.0989400395497483, -0.0145213144685405),
        (36, -0.0979528917494459, -0.0194840415895657),
        (37, -0.0969657439491436, -0.0244467687105909),
        (38, -0.0914086774858697, -0.0378627139332354),
        (39, -0.0858516110225958, -0.05127865915588),
        (40, -0.0830404519257676, -0.0554858560599625),
        (41, -0.0802292928289395, -0.0596930529640449),
        (42, -0.0699611728964922, -0.0699611728964922),
        (43, -0.059693052964045, -0.0802292928289395),
        (44, -0.0554858560599625, -0.0830404519257676),
        (45, -0.05127865915588, -0.0858516110225958),
        (46, -0.0378627139332354, -0.0914086774858697),
        (47, -0.0244467687105908, -0.0969657439491436),
        (48, -0.0194840415895657, -0.0979528917494459),
        (49, -0.0145213144685405, -0.0989400395497483),
        (50, 0.0, -0.0989400395497483),
        (51, 0.0145213144685405, -0.0989400395497483),
        (52, 0.0194840415895657, -0.0979528917494459),
        (53, 0.0244467687105909, -0.0969657439491436),
        (54, 0.0378627139332354, -0.0914086774858697),
        (55, 0.05127865915588, -0.0858516110225958),
        (56, 0.0554858560599625, -0.0830404519257676),
        (57, 0.0596930529640449, -0.0802292928289395),
        (58, 0.0699611728964922, -0.0699611728964922),
        (59, 0.0802292928289396, -0.0596930529640449),
        (60, 0.0830404519257677, -0.0554858560599624),
        (61, 0.0858516110225958, -0.0512786591558799),
        (62, 0.0914086774858697, -0.0378627139332354),
        (63, 0.0969657439491436, -0.0244467687105908),
        (64, 0.0979528917494459, -0.0194840415895657),
    ]
    for tag, x, y in outer_coords:
        p(x, y, lc, tag)
    for i in range(1, 65):
        gmsh.model.geo.addLine(i, i + 1 if i < 64 else 1, 100 + i)

    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 1000)

    # Center
    p(0, 0, lc2, 2000)

    rid = {"r0": r0, "r1": r1, "r2": r2, "r3": r3}
    circle_loops = {}
    offset = 3000

    for i, (rk, rad) in enumerate(sorted(rid.items(), key=lambda item: item[1])):
        start = offset + i * 10
        tags = {"E": start + 1, "N": start + 2, "W": start + 3, "S": start + 4}

        p(rad * cth, rad * sth, lc2, tags["E"])
        p(-rad * sth, rad * cth, lc2, tags["N"])
        p(-rad * cth, -rad * sth, lc2, tags["W"])
        p(rad * sth, -rad * cth, lc2, tags["S"])

        arcs = [start + 5, start + 6, start + 7, start + 8]
        circle(tags["E"], 2000, tags["N"], arcs[0])
        circle(tags["N"], 2000, tags["W"], arcs[1])
        circle(tags["W"], 2000, tags["S"], arcs[2])
        circle(tags["S"], 2000, tags["E"], arcs[3])

        loop_id = 1100 + i
        gmsh.model.geo.addCurveLoop(arcs, loop_id)
        circle_loops[rk] = loop_id

    boundary_lines = {
        "PEE1": [101, 102],
        "PEE2": [162, 161],
        "PEE3": [158, 157],
        "PEE4": [154, 153],
        "PEE5": [150, 149],
        "PEE6": [146, 145],
        "PEE7": [142, 141],
        "PEE8": [138, 137],
        "PEE1a": [134, 133],
        "PEE2a": [130, 129],
        "PEE3a": [126, 125],
        "PEE4a": [122, 121],
        "PEE5a": [118, 117],
        "PEE6a": [114, 113],
        "PEE7a": [110, 109],
        "PEE8a": [106, 105],
        "DSS1": [163, 164],
        "DSS2": [160, 159],
        "DSS3": [156, 155],
        "DSS4": [152, 151],
        "DSS5": [148, 147],
        "DSS6": [144, 143],
        "DSS7": [140, 139],
        "DSS8": [136, 135],
        "DSS1a": [132, 131],
        "DSS2a": [128, 127],
        "DSS3a": [124, 123],
        "DSS4a": [120, 119],
        "DSS5a": [116, 115],
        "DSS6a": [112, 111],
        "DSS7a": [108, 107],
        "DSS8a": [104, 103],
    }

    surface_id = 200
    radii_sorted = sorted(rid, key=rid.get)
    gmsh.model.geo.addPlaneSurface([circle_loops[radii_sorted[0]]], surface_id)

    for i in range(1, len(radii_sorted)):
        surface_id += 1
        gmsh.model.geo.addPlaneSurface(
            [circle_loops[radii_sorted[i]], circle_loops[radii_sorted[i - 1]]],
            surface_id,
        )

    # Outer boundary
    surface_id += 1
    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 3000)
    gmsh.model.geo.addPlaneSurface([3000, circle_loops[radii_sorted[-1]]], surface_id)

    # Recombine surfaces
    for sid in range(200, 204):
        gmsh.model.geo.mesh.setRecombine(2, sid)

    p(0.04, 0.0, lc2, 700)
    gmsh.model.geo.addLine(2000, 700, 1234)  # Fracture line

    gmsh.model.geo.synchronize()

    # Physical Embed fracture line
    gmsh.model.mesh.embed(1, [1234], 2, 200)
    # Physical Source point at center
    gmsh.model.mesh.embed(0, [2000], 2, 200)
    gmsh.model.mesh.embed(0, [700], 2, 200)
    gmsh.model.mesh.generate(2)

    if mode == "BC":
        for name, lines in boundary_lines.items():
            gmsh.model.addPhysicalGroup(1, lines, name=name)
        gmsh.model.addPhysicalGroup(0, [2000], name="Source")
        gmsh.model.addPhysicalGroup(0, [700], name="Outlet")
        gmsh.model.addPhysicalGroup(0, [2], name="p_right")
        gmsh.model.addPhysicalGroup(0, [18], name="p_top")
        gmsh.model.addPhysicalGroup(0, [34], name="p_left")
        gmsh.model.addPhysicalGroup(0, [50], name="p_bottom")
    elif mode == "domain":
        gmsh.model.addPhysicalGroup(2, [200], name="Central_sample")
        gmsh.model.addPhysicalGroup(2, [201, 202, 203], name="OuterPart_sample")
        gmsh.model.addPhysicalGroup(2, [204], name="Rubber_sheath")
        gmsh.model.addPhysicalGroup(1, [1234], name="Fracture")

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    msh_file = Path(out_dir, f"{meshname}.msh")
    gmsh.write(str(msh_file))
    gmsh.finalize()
    return msh_file


def mesh_GreatCell_fullFracture(
    lc=0.0025,
    lc2=0.0025,
    r0=0.097,
    r1=0.094,
    r2=0.09,
    r3=0.065,
    fracture_angle_deg=0.0,
    out_dir=".",
    meshname="Greatcell_mesh",
    mode="domain",
):
    out_path = Path(out_dir)
    if not out_path.exists():
        out_path.mkdir(parents=True)
        print(f"Created output directory: {out_dir}")
    else:
        print(f"Output directory already exists: {out_dir}")

    if not gmsh.isInitialized():
        gmsh.initialize()

    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.model.add("circular_fracture")

    th = math.radians(fracture_angle_deg)
    sth = math.sin(th)
    cth = math.cos(th)

    def p(x, y, mesh_size, tag):
        gmsh.model.geo.addPoint(x, y, 0, mesh_size, tag)

    def circle(p1, pc, p2, tag):
        gmsh.model.geo.addCircleArc(p1, pc, p2, tag)

    # Outer points
    outer_coords = [
        (1, 0.0989400395497483, -0.0145213144685405),
        (2, 0.0989400395497483, 0.0),
        (3, 0.0989400395497483, 0.0145213144685405),
        (4, 0.0979528917494459, 0.0194840415895657),
        (5, 0.0969657439491436, 0.0244467687105909),
        (6, 0.0914086774858697, 0.0378627139332354),
        (7, 0.0858516110225958, 0.05127865915588),
        (8, 0.0830404519257676, 0.0554858560599625),
        (9, 0.0802292928289395, 0.0596930529640449),
        (10, 0.0699611728964922, 0.0699611728964922),
        (11, 0.0596930529640449, 0.0802292928289395),
        (12, 0.0554858560599625, 0.0830404519257676),
        (13, 0.05127865915588, 0.0858516110225958),
        (14, 0.0378627139332354, 0.0914086774858697),
        (15, 0.0244467687105908, 0.0969657439491436),
        (16, 0.0194840415895657, 0.0979528917494459),
        (17, 0.0145213144685405, 0.0989400395497483),
        (18, 0.0, 0.0989400395497483),
        (19, -0.0145213144685405, 0.0989400395497483),
        (20, -0.0194840415895657, 0.0979528917494459),
        (21, -0.0244467687105908, 0.0969657439491436),
        (22, -0.0378627139332354, 0.0914086774858697),
        (23, -0.05127865915588, 0.0858516110225958),
        (24, -0.0554858560599625, 0.0830404519257676),
        (25, -0.0596930529640449, 0.0802292928289395),
        (26, -0.0699611728964922, 0.0699611728964922),
        (27, -0.0802292928289395, 0.0596930529640449),
        (28, -0.0830404519257677, 0.0554858560599624),
        (29, -0.0858516110225958, 0.05127865915588),
        (30, -0.0914086774858697, 0.0378627139332354),
        (31, -0.0969657439491436, 0.0244467687105908),
        (32, -0.0979528917494459, 0.0194840415895657),
        (33, -0.0989400395497483, 0.0145213144685405),
        (34, -0.0989400395497483, 0.0),
        (35, -0.0989400395497483, -0.0145213144685405),
        (36, -0.0979528917494459, -0.0194840415895657),
        (37, -0.0969657439491436, -0.0244467687105909),
        (38, -0.0914086774858697, -0.0378627139332354),
        (39, -0.0858516110225958, -0.05127865915588),
        (40, -0.0830404519257676, -0.0554858560599625),
        (41, -0.0802292928289395, -0.0596930529640449),
        (42, -0.0699611728964922, -0.0699611728964922),
        (43, -0.059693052964045, -0.0802292928289395),
        (44, -0.0554858560599625, -0.0830404519257676),
        (45, -0.05127865915588, -0.0858516110225958),
        (46, -0.0378627139332354, -0.0914086774858697),
        (47, -0.0244467687105908, -0.0969657439491436),
        (48, -0.0194840415895657, -0.0979528917494459),
        (49, -0.0145213144685405, -0.0989400395497483),
        (50, 0.0, -0.0989400395497483),
        (51, 0.0145213144685405, -0.0989400395497483),
        (52, 0.0194840415895657, -0.0979528917494459),
        (53, 0.0244467687105909, -0.0969657439491436),
        (54, 0.0378627139332354, -0.0914086774858697),
        (55, 0.05127865915588, -0.0858516110225958),
        (56, 0.0554858560599625, -0.0830404519257676),
        (57, 0.0596930529640449, -0.0802292928289395),
        (58, 0.0699611728964922, -0.0699611728964922),
        (59, 0.0802292928289396, -0.0596930529640449),
        (60, 0.0830404519257677, -0.0554858560599624),
        (61, 0.0858516110225958, -0.0512786591558799),
        (62, 0.0914086774858697, -0.0378627139332354),
        (63, 0.0969657439491436, -0.0244467687105908),
        (64, 0.0979528917494459, -0.0194840415895657),
    ]
    for tag, x, y in outer_coords:
        p(x, y, lc, tag)
    for i in range(1, 65):
        gmsh.model.geo.addLine(i, i + 1 if i < 64 else 1, 100 + i)

    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 1000)
    # Center point
    p(0, 0, lc2, 2000)

    # r3 (smallest)
    p(r3 * cth, r3 * sth, lc2, 3001)
    p(-r3 * sth, r3 * cth, lc2, 3002)
    p(-r3 * cth, -r3 * sth, lc2, 3003)
    p(r3 * sth, -r3 * cth, lc2, 3004)
    circle(3001, 2000, 3002, 3005)
    circle(3002, 2000, 3003, 3006)
    circle(3003, 2000, 3004, 3007)
    circle(3004, 2000, 3001, 3008)

    # r2
    p(r2 * cth, r2 * sth, lc2, 3011)
    p(-r2 * sth, r2 * cth, lc2, 3012)
    p(-r2 * cth, -r2 * sth, lc2, 3013)
    p(r2 * sth, -r2 * cth, lc2, 3014)
    circle(3011, 2000, 3012, 3015)
    circle(3012, 2000, 3013, 3016)
    circle(3013, 2000, 3014, 3017)
    circle(3014, 2000, 3011, 3018)

    # r1
    p(r1 * cth, r1 * sth, lc2, 3021)
    p(-r1 * sth, r1 * cth, lc2, 3022)
    p(-r1 * cth, -r1 * sth, lc2, 3023)
    p(r1 * sth, -r1 * cth, lc2, 3024)
    circle(3021, 2000, 3022, 3025)
    circle(3022, 2000, 3023, 3026)
    circle(3023, 2000, 3024, 3027)
    circle(3024, 2000, 3021, 3028)

    # r0 (largest inner circle)
    p(r0 * cth, r0 * sth, lc2, 3031)
    p(-r0 * sth, r0 * cth, lc2, 3032)
    p(-r0 * cth, -r0 * sth, lc2, 3033)
    p(r0 * sth, -r0 * cth, lc2, 3034)
    circle(3031, 2000, 3032, 3035)
    circle(3032, 2000, 3033, 3036)
    circle(3033, 2000, 3034, 3037)
    circle(3034, 2000, 3031, 3038)

    split_lines = []

    gmsh.model.geo.addLine(2000, 3001, 8000)  # center to r3-right
    split_lines.append(8000)

    gmsh.model.geo.addLine(3001, 3011, 8001)  # r3-right to r2-right
    split_lines.append(8001)

    gmsh.model.geo.addLine(3011, 3021, 8002)  # r2-right to r1-right
    split_lines.append(8002)

    gmsh.model.geo.addLine(3021, 3031, 8003)  # r1-right to r0-right
    split_lines.append(8003)

    gmsh.model.geo.addLine(3033, 3023, 8004)  # r0-left to r1-left
    split_lines.append(8004)

    gmsh.model.geo.addLine(3023, 3013, 8005)  # r1-left to r2-left
    split_lines.append(8005)

    gmsh.model.geo.addLine(3013, 3003, 8006)  # r2-left to r3-left
    split_lines.append(8006)

    gmsh.model.geo.addLine(3003, 2000, 8007)  # r3-left to center
    split_lines.append(8007)

    # Top halves
    gmsh.model.geo.addCurveLoop([3005, 3006, 8007, 8000], 200)
    gmsh.model.geo.addPlaneSurface([200], 200)

    gmsh.model.geo.addCurveLoop([3015, 3016, 8006, -3006, -3005, 8001], 201)
    gmsh.model.geo.addPlaneSurface([201], 201)

    gmsh.model.geo.addCurveLoop([3025, 3026, 8005, -3016, -3015, 8002], 202)
    gmsh.model.geo.addPlaneSurface([202], 202)

    gmsh.model.geo.addCurveLoop([3035, 3036, 8004, -3026, -3025, 8003], 203)
    gmsh.model.geo.addPlaneSurface([203], 203)

    # Bottom halves
    gmsh.model.geo.addCurveLoop([3007, 3008, -8000, -8007], 204)
    gmsh.model.geo.addPlaneSurface([204], 204)

    gmsh.model.geo.addCurveLoop([3017, 3018, -8001, -3008, -3007, -8006], 205)
    gmsh.model.geo.addPlaneSurface([205], 205)

    gmsh.model.geo.addCurveLoop([3027, 3028, -8002, -3018, -3017, -8005], 206)
    gmsh.model.geo.addPlaneSurface([206], 206)

    gmsh.model.geo.addCurveLoop([3037, 3038, -8003, -3028, -3027, -8004], 207)
    gmsh.model.geo.addPlaneSurface([207], 207)

    # Outer rubber sheath (between outer boundary and r0 circle)
    gmsh.model.geo.addCurveLoop([3035, 3036, 3037, 3038], 208)
    gmsh.model.geo.addPlaneSurface([1000, 208], 208)

    boundary_lines = {
        "PEE1": [101, 102],
        "PEE2": [162, 161],
        "PEE3": [158, 157],
        "PEE4": [154, 153],
        "PEE5": [150, 149],
        "PEE6": [146, 145],
        "PEE7": [142, 141],
        "PEE8": [138, 137],
        "PEE1a": [134, 133],
        "PEE2a": [130, 129],
        "PEE3a": [126, 125],
        "PEE4a": [122, 121],
        "PEE5a": [118, 117],
        "PEE6a": [114, 113],
        "PEE7a": [110, 109],
        "PEE8a": [106, 105],
        "DSS1": [163, 164],
        "DSS2": [160, 159],
        "DSS3": [156, 155],
        "DSS4": [152, 151],
        "DSS5": [148, 147],
        "DSS6": [144, 143],
        "DSS7": [140, 139],
        "DSS8": [136, 135],
        "DSS1a": [132, 131],
        "DSS2a": [128, 127],
        "DSS3a": [124, 123],
        "DSS4a": [120, 119],
        "DSS5a": [116, 115],
        "DSS6a": [112, 111],
        "DSS7a": [108, 107],
        "DSS8a": [104, 103],
    }

    # Recombine surfaces
    for sid in range(200, 208):
        gmsh.model.geo.mesh.setRecombine(2, sid)

    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    if mode == "BC":
        for name, lines in boundary_lines.items():
            gmsh.model.addPhysicalGroup(1, lines, name=name)
        gmsh.model.addPhysicalGroup(0, [2000], name="Source")
        gmsh.model.addPhysicalGroup(0, [2], name="p_right")
        gmsh.model.addPhysicalGroup(0, [18], name="p_top")
        gmsh.model.addPhysicalGroup(0, [34], name="p_left")
        gmsh.model.addPhysicalGroup(0, [50], name="p_bottom")
        gmsh.model.addPhysicalGroup(0, [3031], name="outlet1")
        gmsh.model.addPhysicalGroup(0, [3032], name="outlet2")

    elif mode == "domain":
        # Physical surfaces
        gmsh.model.addPhysicalGroup(2, [200, 204], name="Central_sample")
        gmsh.model.addPhysicalGroup(
            2, [201, 202, 203, 205, 206, 207], name="OuterPart_sample"
        )
        gmsh.model.addPhysicalGroup(2, [208], name="Rubber_sheath")
        gmsh.model.addPhysicalGroup(
            1, [8000, 8001, 8002, 8003, 8004, 8005, 8006, 8007], name="fracture"
        )

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    msh_file = Path(out_dir, f"{meshname}.msh")
    gmsh.write(str(msh_file))
    gmsh.finalize()
    return msh_file


def mesh_GreatCell_VPF(
    lc2=0.0005,
    r0=0.097,
    r1=0.090,
    r2=0.065,
    r3=0.0,
    out_dir=".",
    meshname="GreatCell_mesh",
    mode="domain",
):
    out_path = Path(out_dir)
    if not out_path.exists():
        out_path.mkdir(parents=True)
        print(f"Created output directory: {out_dir}")
    else:
        print(f"Output directory already exists: {out_dir}")

    if not gmsh.isInitialized():
        gmsh.initialize()

    if not gmsh.isInitialized():
        gmsh.initialize()
    else:
        gmsh.clear()

    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.model.add("circular_fracture")

    crack_thickness = 0.5 * lc2
    fracture_angle_to_bag1_sx = 0.0
    frac = 0.04

    th = fracture_angle_to_bag1_sx * math.pi / 180
    sth = math.sin(th)
    cth = math.cos(th)

    def p(x, y, mesh_size, tag):
        gmsh.model.geo.addPoint(x, y, 0, mesh_size, tag)

    def line(p1, p2, tag):
        gmsh.model.geo.addLine(p1, p2, tag)

    def circle(p1, pc, p2, tag):
        gmsh.model.geo.addCircleArc(p1, pc, p2, tag)

    # Outer points
    outer_coords = [
        (1, 0.0989400395497483, -0.0145213144685405),
        (2, 0.0989400395497483, 0.0),
        (3, 0.0989400395497483, 0.0145213144685405),
        (4, 0.0979528917494459, 0.0194840415895657),
        (5, 0.0969657439491436, 0.0244467687105909),
        (6, 0.0914086774858697, 0.0378627139332354),
        (7, 0.0858516110225958, 0.05127865915588),
        (8, 0.0830404519257676, 0.0554858560599625),
        (9, 0.0802292928289395, 0.0596930529640449),
        (10, 0.0699611728964922, 0.0699611728964922),
        (11, 0.0596930529640449, 0.0802292928289395),
        (12, 0.0554858560599625, 0.0830404519257676),
        (13, 0.05127865915588, 0.0858516110225958),
        (14, 0.0378627139332354, 0.0914086774858697),
        (15, 0.0244467687105908, 0.0969657439491436),
        (16, 0.0194840415895657, 0.0979528917494459),
        (17, 0.0145213144685405, 0.0989400395497483),
        (18, 0.0, 0.0989400395497483),
        (19, -0.0145213144685405, 0.0989400395497483),
        (20, -0.0194840415895657, 0.0979528917494459),
        (21, -0.0244467687105908, 0.0969657439491436),
        (22, -0.0378627139332354, 0.0914086774858697),
        (23, -0.05127865915588, 0.0858516110225958),
        (24, -0.0554858560599625, 0.0830404519257676),
        (25, -0.0596930529640449, 0.0802292928289395),
        (26, -0.0699611728964922, 0.0699611728964922),
        (27, -0.0802292928289395, 0.0596930529640449),
        (28, -0.0830404519257677, 0.0554858560599624),
        (29, -0.0858516110225958, 0.05127865915588),
        (30, -0.0914086774858697, 0.0378627139332354),
        (31, -0.0969657439491436, 0.0244467687105908),
        (32, -0.0979528917494459, 0.0194840415895657),
        (33, -0.0989400395497483, 0.0145213144685405),
        (34, -0.0989400395497483, 0.0),
        (35, -0.0989400395497483, -0.0145213144685405),
        (36, -0.0979528917494459, -0.0194840415895657),
        (37, -0.0969657439491436, -0.0244467687105909),
        (38, -0.0914086774858697, -0.0378627139332354),
        (39, -0.0858516110225958, -0.05127865915588),
        (40, -0.0830404519257676, -0.0554858560599625),
        (41, -0.0802292928289395, -0.0596930529640449),
        (42, -0.0699611728964922, -0.0699611728964922),
        (43, -0.059693052964045, -0.0802292928289395),
        (44, -0.0554858560599625, -0.0830404519257676),
        (45, -0.05127865915588, -0.0858516110225958),
        (46, -0.0378627139332354, -0.0914086774858697),
        (47, -0.0244467687105908, -0.0969657439491436),
        (48, -0.0194840415895657, -0.0979528917494459),
        (49, -0.0145213144685405, -0.0989400395497483),
        (50, 0.0, -0.0989400395497483),
        (51, 0.0145213144685405, -0.0989400395497483),
        (52, 0.0194840415895657, -0.0979528917494459),
        (53, 0.0244467687105909, -0.0969657439491436),
        (54, 0.0378627139332354, -0.0914086774858697),
        (55, 0.05127865915588, -0.0858516110225958),
        (56, 0.0554858560599625, -0.0830404519257676),
        (57, 0.0596930529640449, -0.0802292928289395),
        (58, 0.0699611728964922, -0.0699611728964922),
        (59, 0.0802292928289396, -0.0596930529640449),
        (60, 0.0830404519257677, -0.0554858560599624),
        (61, 0.0858516110225958, -0.0512786591558799),
        (62, 0.0914086774858697, -0.0378627139332354),
        (63, 0.0969657439491436, -0.0244467687105908),
        (64, 0.0979528917494459, -0.0194840415895657),
    ]
    for tag, x, y in outer_coords:
        p(x, y, lc2, tag)

    for i in range(1, 65):
        line(i, i + 1 if i < 64 else 1, 100 + i)

    gmsh.model.geo.addCurveLoop(list(range(101, 165)), 1000)

    p(0.0, 0.0, lc2, 2000)

    # Define fracture-related points (r0, r1, r2, etc.)
    p3011 = gmsh.model.geo.addPoint(
        frac * cth - crack_thickness * sth,
        frac * sth + crack_thickness * cth,
        0,
        lc2,
        3011,
    )
    p3012 = gmsh.model.geo.addPoint(
        frac * cth + crack_thickness * sth,
        frac * sth - crack_thickness * cth,
        0,
        lc2,
        3012,
    )

    p3001 = gmsh.model.geo.addPoint(
        r0 * cth - crack_thickness * sth, r0 * sth + crack_thickness * cth, 0, lc2, 3001
    )
    p3002 = gmsh.model.geo.addPoint(
        r0 * cth + crack_thickness * sth, r0 * sth - crack_thickness * cth, 0, lc2, 3002
    )
    p3003 = gmsh.model.geo.addPoint(
        -r0 * cth - crack_thickness * sth,
        -r0 * sth + crack_thickness * cth,
        0,
        lc2,
        3003,
    )
    p3004 = gmsh.model.geo.addPoint(
        -r0 * cth + crack_thickness * sth,
        -r0 * sth - crack_thickness * cth,
        0,
        lc2,
        3004,
    )
    p3005 = gmsh.model.geo.addPoint(-r0 * sth, r0 * cth, 0, lc2, 3005)
    p3006 = gmsh.model.geo.addPoint(r0 * sth, -r0 * cth, 0, lc2, 3006)

    p3101 = gmsh.model.geo.addPoint(
        r1 * cth - crack_thickness * sth, r1 * sth + crack_thickness * cth, 0, lc2, 3101
    )
    p3102 = gmsh.model.geo.addPoint(
        r1 * cth + crack_thickness * sth, r1 * sth - crack_thickness * cth, 0, lc2, 3102
    )
    p3103 = gmsh.model.geo.addPoint(
        -r1 * cth - crack_thickness * sth,
        -r1 * sth + crack_thickness * cth,
        0,
        lc2,
        3103,
    )
    p3104 = gmsh.model.geo.addPoint(
        -r1 * cth + crack_thickness * sth,
        -r1 * sth - crack_thickness * cth,
        0,
        lc2,
        3104,
    )
    p3105 = gmsh.model.geo.addPoint(-r1 * sth, r1 * cth, 0, lc2, 3105)
    p3106 = gmsh.model.geo.addPoint(r1 * sth, -r1 * cth, 0, lc2, 3106)

    p3201 = gmsh.model.geo.addPoint(
        r2 * cth - crack_thickness * sth, r2 * sth + crack_thickness * cth, 0, lc2, 3201
    )
    p3202 = gmsh.model.geo.addPoint(
        r2 * cth + crack_thickness * sth, r2 * sth - crack_thickness * cth, 0, lc2, 3202
    )
    p3203 = gmsh.model.geo.addPoint(
        -r2 * cth - crack_thickness * sth,
        -r2 * sth + crack_thickness * cth,
        0,
        lc2,
        3203,
    )
    p3204 = gmsh.model.geo.addPoint(
        -r2 * cth + crack_thickness * sth,
        -r2 * sth - crack_thickness * cth,
        0,
        lc2,
        3204,
    )
    p3205 = gmsh.model.geo.addPoint(-r2 * sth, r2 * cth, 0, lc2, 3205)
    p3206 = gmsh.model.geo.addPoint(r2 * sth, -r2 * cth, 0, lc2, 3206)

    p3301 = gmsh.model.geo.addPoint(
        r3 * cth - crack_thickness * sth, r3 * sth + crack_thickness * cth, 0, lc2, 3301
    )
    p3302 = gmsh.model.geo.addPoint(
        r3 * cth + crack_thickness * sth, r3 * sth - crack_thickness * cth, 0, lc2, 3302
    )
    p3303 = gmsh.model.geo.addPoint(
        -r3 * cth - crack_thickness * sth,
        -r3 * sth + crack_thickness * cth,
        0,
        lc2,
        3303,
    )
    p3304 = gmsh.model.geo.addPoint(
        -r3 * cth + crack_thickness * sth,
        -r3 * sth - crack_thickness * cth,
        0,
        lc2,
        3304,
    )
    p3305 = gmsh.model.geo.addPoint(-r3 * sth, r3 * cth, 0, lc2, 3305)
    p3306 = gmsh.model.geo.addPoint(r3 * sth, -r3 * cth, 0, lc2, 3306)

    p3501 = gmsh.model.geo.addPoint(
        -crack_thickness * sth, crack_thickness * cth, 0, lc2, 3501
    )
    p3502 = gmsh.model.geo.addPoint(
        crack_thickness * sth, -crack_thickness * cth, 0, lc2, 3502
    )

    # # Arcs and lines for fracture
    c401 = circle(p3001, 2000, p3005, 401)
    c402 = circle(p3005, 2000, p3003, 402)
    c403 = circle(p3004, 2000, p3006, 403)
    c404 = circle(p3006, 2000, p3002, 404)

    c411 = circle(p3101, 2000, p3105, 411)
    c412 = circle(p3105, 2000, p3103, 412)
    c413 = circle(p3104, 2000, p3106, 413)
    c414 = circle(p3106, 2000, p3102, 414)

    c421 = circle(p3201, 2000, p3205, 421)
    c422 = circle(p3205, 2000, p3203, 422)
    c423 = circle(p3204, 2000, p3206, 423)
    c424 = circle(p3206, 2000, p3202, 424)

    c431 = circle(p3301, 2000, p3305, 431)
    c432 = circle(p3305, 2000, p3303, 432)
    c433 = circle(p3304, 2000, p3306, 433)
    c434 = circle(p3306, 2000, p3302, 434)

    # Lines for fracture
    l501 = line(p3001, p3101, 501)
    l502 = line(p3101, p3201, 502)
    l503 = line(p3201, p3301, 503)

    l504 = line(p3301, p3011, 504)
    l5040 = line(p3011, p3501, 5040)
    l505 = line(p3501, p3303, 505)

    l506 = line(p3303, p3203, 506)
    l507 = line(p3203, p3103, 507)
    l508 = line(p3103, p3003, 508)

    l509 = line(p3002, p3102, 509)
    l510 = line(p3102, p3202, 510)
    l511 = line(p3202, p3302, 511)

    l512 = line(p3302, p3012, 512)
    l5120 = line(p3012, p3502, 5120)

    l513 = line(p3502, p3304, 513)

    l514 = line(p3304, p3204, 514)
    l515 = line(p3204, p3104, 515)
    l516 = line(p3104, p3004, 516)

    l517 = line(p3003, p3004, 517)
    l518 = line(p3103, p3104, 518)
    l519 = line(p3203, p3204, 519)
    l520 = line(p3303, p3304, 520)

    l521 = line(p3501, p3502, 521)
    l526 = line(p3011, p3012, 526)

    l522 = line(p3001, p3002, 522)
    l523 = line(p3101, p3102, 523)
    l524 = line(p3201, p3202, 524)
    l525 = line(p3301, p3302, 525)

    # Define Curve Loops
    gmsh.model.geo.addCurveLoop([-l525, l504, l526, -l512], 601)
    gmsh.model.geo.addCurveLoop([l5040, l521, -l5120, -l526], 6010)

    gmsh.model.geo.addCurveLoop([l521, l513, -l520, -l505], 602)

    gmsh.model.geo.addCurveLoop([l504, l5040, l505, -c432, -c431], 603)
    gmsh.model.geo.addCurveLoop([l512, l5120, l513, c433, c434], 604)

    gmsh.model.geo.addCurveLoop([l506, l519, -l514, -l520], 605)
    gmsh.model.geo.addCurveLoop([l503, l525, -l511, -l524], 606)

    gmsh.model.geo.addCurveLoop([c432, c431, l506, -c422, -c421, l503], 607)
    gmsh.model.geo.addCurveLoop([-c434, -c433, l514, c423, c424, l511], 608)

    gmsh.model.geo.addCurveLoop([l507, l518, -l515, -l519], 609)
    gmsh.model.geo.addCurveLoop([l502, l524, -l510, -l523], 610)

    gmsh.model.geo.addCurveLoop([c422, c421, l507, -c412, -c411, l502], 611)
    gmsh.model.geo.addCurveLoop([-c424, -c423, l515, c413, c414, l510], 612)

    gmsh.model.geo.addCurveLoop([l508, l517, -l516, -l518], 613)
    gmsh.model.geo.addCurveLoop([l501, l523, -l509, -l522], 614)

    gmsh.model.geo.addCurveLoop([c412, c411, l508, -c402, -c401, l501], 615)
    gmsh.model.geo.addCurveLoop([-c414, -c413, l516, c403, c404, l509], 616)

    gmsh.model.geo.addCurveLoop([c401, c402, l517, c403, c404, -l522], 617)

    # # Plane Surfaces
    gmsh.model.geo.addPlaneSurface([601], 701)
    gmsh.model.geo.addPlaneSurface([6010], 702)
    gmsh.model.geo.addPlaneSurface([602], 703)
    gmsh.model.geo.addPlaneSurface([603], 704)
    gmsh.model.geo.addPlaneSurface([604], 705)
    gmsh.model.geo.addPlaneSurface([605], 706)
    gmsh.model.geo.addPlaneSurface([606], 707)
    gmsh.model.geo.addPlaneSurface([607], 708)
    gmsh.model.geo.addPlaneSurface([608], 709)
    gmsh.model.geo.addPlaneSurface([609], 710)
    gmsh.model.geo.addPlaneSurface([610], 711)
    gmsh.model.geo.addPlaneSurface([611], 712)
    gmsh.model.geo.addPlaneSurface([612], 713)
    gmsh.model.geo.addPlaneSurface([613], 714)
    gmsh.model.geo.addPlaneSurface([614], 715)
    gmsh.model.geo.addPlaneSurface([615], 716)
    gmsh.model.geo.addPlaneSurface([616], 717)
    gmsh.model.geo.addPlaneSurface([617, 1000], 718)

    # Transfinite settings
    for curve_id in [501, 509, 508, 516]:
        gmsh.model.geo.mesh.setTransfiniteCurve(
            curve_id, max(1, int((r0 - r1) / lc2)) + 1
        )
    for curve_id in [502, 510, 515, 507]:
        gmsh.model.geo.mesh.setTransfiniteCurve(
            curve_id, max(1, int((r1 - r2) / lc2)) + 1
        )
    for curve_id in [503, 511, 514, 506]:
        gmsh.model.geo.mesh.setTransfiniteCurve(
            curve_id, max(1, int((r2 - r3) / lc2)) + 1
        )

    for curve_id in [504, 512]:
        gmsh.model.geo.mesh.setTransfiniteCurve(curve_id, int((r3 - frac) / lc2) + 1)

    for curve_id in [505, 513]:
        gmsh.model.geo.mesh.setTransfiniteCurve(curve_id, int(r3 / lc2))

    for curve_id in [5040, 5120]:
        gmsh.model.geo.mesh.setTransfiniteCurve(curve_id, int(frac / lc2))

    for curve_id in [517, 518, 519, 520, 521, 522, 523, 524, 525, 526]:
        gmsh.model.geo.mesh.setTransfiniteCurve(curve_id, 1)

    boundary_lines = {
        "PEE1": [101, 102],
        "PEE2": [162, 161],
        "PEE3": [158, 157],
        "PEE4": [154, 153],
        "PEE5": [150, 149],
        "PEE6": [146, 145],
        "PEE7": [142, 141],
        "PEE8": [138, 137],
        "PEE1a": [134, 133],
        "PEE2a": [130, 129],
        "PEE3a": [126, 125],
        "PEE4a": [122, 121],
        "PEE5a": [118, 117],
        "PEE6a": [114, 113],
        "PEE7a": [110, 109],
        "PEE8a": [106, 105],
        "DSS1": [163, 164],
        "DSS2": [160, 159],
        "DSS3": [156, 155],
        "DSS4": [152, 151],
        "DSS5": [148, 147],
        "DSS6": [144, 143],
        "DSS7": [140, 139],
        "DSS8": [136, 135],
        "DSS1a": [132, 131],
        "DSS2a": [128, 127],
        "DSS3a": [124, 123],
        "DSS4a": [120, 119],
        "DSS5a": [116, 115],
        "DSS6a": [112, 111],
        "DSS7a": [108, 107],
        "DSS8a": [104, 103],
    }

    # Surface recombination for square elements
    for sid in range(700, 719):
        gmsh.model.geo.mesh.setRecombine(2, sid)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)

    if mode == "BC":
        for name, lines in boundary_lines.items():
            gmsh.model.addPhysicalGroup(1, lines, name=name)
        gmsh.model.addPhysicalGroup(0, [2], name="p_right")
        gmsh.model.addPhysicalGroup(0, [18], name="p_top")
        gmsh.model.addPhysicalGroup(0, [34], name="p_left")
        gmsh.model.addPhysicalGroup(0, [50], name="p_bottom")
        gmsh.model.addPhysicalGroup(1, [521], name="SourceLine")
        gmsh.model.addPhysicalGroup(1, [526], name="OutletLine_embeddedFracture")
        gmsh.model.addPhysicalGroup(1, [523], name="OutletLine_fullFracture_right")
        gmsh.model.addPhysicalGroup(1, [518], name="OutletLine_fullFracture_left")

    elif mode == "domain":
        # Physical surfaces
        gmsh.model.addPhysicalGroup(
            2,
            [
                701,
                702,
                703,
                704,
                705,
                706,
                707,
                708,
                709,
                710,
                711,
                712,
                713,
                714,
                715,
                716,
                717,
            ],
            name="Central_sample",
        )
        gmsh.model.addPhysicalGroup(2, [718], name="Rubber_sheath")

    msh_file = Path(out_dir, f"{meshname}.msh")
    gmsh.write(str(msh_file))
    gmsh.finalize()
    return msh_file
