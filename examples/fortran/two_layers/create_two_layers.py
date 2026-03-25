import gmsh
import numpy as np


def build_two_layer_conformal(
    mesh_size=100.0,
    x0=0.0,
    width=10000.0,
    z0=0.0,
    depth=5000.0,
    z_interface=3000.0,
    src=(5000.0, 1500.0),
):
    """
    Conformal 2-layer mesh: interface is a 1D line in the geometry
    (no triangles crossing it).

    Returns:
      xyz      (N,3) points (x,y,z) with y=0, z downward
      tri_conn (M,3) triangles as 0-based point indices
      tri_vel  (M,)  velocity per triangle (CELL_DATA)
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("rect_2layer")

    x1 = x0 + width
    z1 = z0 + depth

    if not (z0 < z_interface < z1):
        raise ValueError("z_interface must be strictly inside the rectangle.")

    sx, sz = src
    if not (x0 <= sx <= x1 and z0 <= sz <= z1):
        raise ValueError("src point must be inside the rectangle.")

    p_bl = gmsh.model.geo.addPoint(x0, z0, 0, mesh_size)
    p_br = gmsh.model.geo.addPoint(x1, z0, 0, mesh_size)
    p_tr = gmsh.model.geo.addPoint(x1, z1, 0, mesh_size)
    p_tl = gmsh.model.geo.addPoint(x0, z1, 0, mesh_size)

    p_il = gmsh.model.geo.addPoint(x0, z_interface, 0, mesh_size)
    p_ir = gmsh.model.geo.addPoint(x1, z_interface, 0, mesh_size)

    p_src = gmsh.model.geo.addPoint(sx, sz, 0, mesh_size)

    l_bottom = gmsh.model.geo.addLine(p_bl, p_br)
    l_right1 = gmsh.model.geo.addLine(p_br, p_ir)
    l_right2 = gmsh.model.geo.addLine(p_ir, p_tr)
    l_top = gmsh.model.geo.addLine(p_tr, p_tl)
    l_left1 = gmsh.model.geo.addLine(p_tl, p_il)
    l_left2 = gmsh.model.geo.addLine(p_il, p_bl)
    l_interface = gmsh.model.geo.addLine(p_il, p_ir)

    loop_top = gmsh.model.geo.addCurveLoop([l_bottom, l_right1, -l_interface, l_left2])
    s_top = gmsh.model.geo.addPlaneSurface([loop_top])

    loop_bot = gmsh.model.geo.addCurveLoop([l_interface, l_right2, l_top, l_left1])
    s_bot = gmsh.model.geo.addPlaneSurface([loop_bot])

    gmsh.model.geo.synchronize()

    gmsh.model.mesh.embed(0, [p_src], 2, s_top)
    gmsh.model.mesh.embed(0, [p_src], 2, s_bot)

    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    xyz = node_coords.reshape(-1, 3).copy()

    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    tri_type = 2
    tri_elem_tags = None
    tri_nodes = None
    for et, tags_block, nodes_block in zip(elem_types, elem_tags, elem_node_tags):
        if et == tri_type:
            tri_elem_tags = np.array(tags_block, dtype=np.int64)
            tri_nodes = np.array(nodes_block, dtype=np.int64).reshape(-1, 3)
            break
    if tri_nodes is None:
        raise RuntimeError("No 3-node triangles found.")

    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}
    tri_conn0 = np.vectorize(lambda t: tag_to_idx[int(t)])(tri_nodes).astype(np.int64)

    tri_vel = np.zeros(len(tri_elem_tags), dtype=float)
    elem_tag_to_row = {int(t): i for i, t in enumerate(tri_elem_tags)}

    def assign_velocity(surface_tag, velocity):
        etypes, etags, _ = gmsh.model.mesh.getElements(2, surface_tag)
        for et, tags_block in zip(etypes, etags):
            if et == tri_type:
                for t in tags_block:
                    tri_vel[elem_tag_to_row[int(t)]] = float(velocity)

    assign_velocity(s_top, 1000.0)
    assign_velocity(s_bot, 3000.0)

    gmsh.finalize()
    return xyz, tri_conn0, tri_vel


def write_vtk_with_point_time(
    filename, xyz, tri_conn0, tri_vel, time_at_points, title="mesh_velocity_time"
):
    npts = xyz.shape[0]
    ntri = tri_conn0.shape[0]
    if len(time_at_points) != npts:
        raise ValueError("time_at_points must have length == number of points")

    total_ints = ntri * 4

    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        f.write(f"POINTS {npts} float\n")
        for p in xyz:
            f.write(f"{p[0]: .8E} {p[1]: .8E} {p[2]: .8E}\n")

        f.write(f"POLYGONS {ntri} {total_ints}\n")
        for t in tri_conn0:
            f.write(f"3 {int(t[0])} {int(t[1])} {int(t[2])}\n")

        f.write(f"\nCELL_DATA {ntri}\n")
        f.write("SCALARS velocity float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for v in tri_vel:
            f.write(f"{float(v):12.5f}\n")

        f.write(f"\nPOINT_DATA {npts}\n")
        f.write("SCALARS time float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for tt in time_at_points:
            f.write(f"{float(tt):12.5f}\n")


def write_vtk_velocity_only(filename, xyz, tri_conn0, tri_vel, title="mesh_velocity"):
    npts = xyz.shape[0]
    ntri = tri_conn0.shape[0]
    total_ints = ntri * 4

    with open(filename, "w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write(f"{title}\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        f.write(f"POINTS {npts} float\n")
        for p in xyz:
            f.write(f"{p[0]: .8E} {p[1]: .8E} {p[2]: .8E}\n")

        f.write(f"POLYGONS {ntri} {total_ints}\n")
        for t in tri_conn0:
            f.write(f"3 {int(t[0])} {int(t[1])} {int(t[2])}\n")

        f.write(f"\nCELL_DATA {ntri}\n")
        f.write("SCALARS velocity float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for v in tri_vel:
            f.write(f"{float(v):12.5f}\n")


def analytical_two_layer_time(xyz, src, z_interface, vup, vdown, dx_search=1.0):
    """
    Analytical first-arrival time for a horizontal two-layer medium.

    Upper-layer nodes use the minimum of:
      - direct arrival in the upper medium,
      - conic/head-wave branch after the critical offset.

    Lower-layer nodes use the refracted branch obtained by minimizing the sum
    of travel times from source to the interface and from the interface to the
    receiver, with the horizontal interface segment between those two offsets.
    """
    sx, sz = src
    ac = np.arcsin(vup / vdown)
    times = np.zeros(xyz.shape[0], dtype=float)
    depth_src = z_interface - sz

    for i, p in enumerate(xyz):
        x = p[0]
        z = p[1]
        x_offset = abs(x - sx)

        if z <= z_interface:
            direct = np.hypot(x_offset, z - sz) / vup
            times[i] = direct

            critical_offset = depth_src * np.tan(ac) + (z_interface - z) * np.tan(ac)
            if x_offset > critical_offset:
                x1 = depth_src * np.tan(ac)
                x2 = (z_interface - z) * np.tan(ac)
                conic = (depth_src + z_interface - z) / (np.cos(ac) * vup)
                conic += (x_offset - x1 - x2) / vdown
                times[i] = min(direct, conic)
        else:
            best = np.hypot(x_offset, z - z_interface) / vdown + depth_src / vup
            offset = 0.0
            while True:
                offset += dx_search
                trial = np.hypot(x_offset - offset, z - z_interface) / vdown
                trial += np.hypot(offset, depth_src) / vup
                if trial < best:
                    best = trial
                else:
                    break
            times[i] = best

    return times


def main():
    mesh_size = 100.0
    x0 = 0.0
    width = 10000.0
    z0 = 0.0
    depth = 5000.0
    z_interface = 3000.0
    vup = 1000.0
    vdown = 3000.0
    src = (5000.0, 1500.0)

    xyz, tri_conn0, tri_vel = build_two_layer_conformal(
        mesh_size=mesh_size,
        x0=x0,
        width=width,
        z0=z0,
        depth=depth,
        z_interface=z_interface,
        src=src,
    )

    theoretical_time = analytical_two_layer_time(
        xyz,
        src=src,
        z_interface=z_interface,
        vup=vup,
        vdown=vdown,
    )

    write_vtk_with_point_time(
        "input_velocity.vtk",
        xyz,
        tri_conn0,
        tri_vel,
        theoretical_time,
        title="Mesh with velocity + analytical point time",
    )

    print("Wrote: input_velocity.vtk")


if __name__ == "__main__":
    main()
