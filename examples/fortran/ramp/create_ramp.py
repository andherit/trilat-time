import gmsh
import numpy as np


# -----------------------------------------------------------------------------
# Geometry / mesh
# -----------------------------------------------------------------------------

def build_complex_two_layer_conformal(
    mesh_size=50.0,
    x0=0.0,
    width=10000.0,
    y0=0.0,
    height=5000.0,
    xk1=3333.3,
    yk1=1666.7,
    xk2=6666.7,
    yk2=2254.42,
    src=(500.0, 500.0),
    vup=1000.0,
    vdown=3000.0,
):
    """
    Build a conformal 2D triangular mesh for a 2-layer model with interface:
      I1: horizontal
      I2: ramp
      I3: horizontal

    Returns
    -------
    xyz : (N, 3) ndarray
        Point coordinates (x, y, z=0).
    tri_conn0 : (M, 3) ndarray
        Triangle connectivity, 0-based.
    tri_vel : (M,) ndarray
        Velocity per triangle.
    """
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("complex_two_layer")

    x1 = x0 + width
    y1 = y0 + height

    sx, sy = src
    if not (x0 <= sx <= x1 and y0 <= sy <= y1):
        raise ValueError("Source must be inside the model domain.")

    # Domain corners
    p_bl = gmsh.model.geo.addPoint(x0, y0, 0.0, mesh_size)
    p_br = gmsh.model.geo.addPoint(x1, y0, 0.0, mesh_size)
    p_tr = gmsh.model.geo.addPoint(x1, y1, 0.0, mesh_size)
    p_tl = gmsh.model.geo.addPoint(x0, y1, 0.0, mesh_size)

    # Interface points
    p_il = gmsh.model.geo.addPoint(x0,  yk1, 0.0, mesh_size)
    p_k1 = gmsh.model.geo.addPoint(xk1, yk1, 0.0, mesh_size)
    p_k2 = gmsh.model.geo.addPoint(xk2, yk2, 0.0, mesh_size)
    p_ir = gmsh.model.geo.addPoint(x1,  yk2, 0.0, mesh_size)

    # Source
    p_src = gmsh.model.geo.addPoint(sx, sy, 0.0, mesh_size)

    # Outer boundary
    l_bottom = gmsh.model.geo.addLine(p_bl, p_br)
    l_right1 = gmsh.model.geo.addLine(p_br, p_ir)
    l_right2 = gmsh.model.geo.addLine(p_ir, p_tr)
    l_top = gmsh.model.geo.addLine(p_tr, p_tl)
    l_left1 = gmsh.model.geo.addLine(p_tl, p_il)
    l_left2 = gmsh.model.geo.addLine(p_il, p_bl)

    # Interface
    l_i1 = gmsh.model.geo.addLine(p_il, p_k1)
    l_i2 = gmsh.model.geo.addLine(p_k1, p_k2)
    l_i3 = gmsh.model.geo.addLine(p_k2, p_ir)

    # Upper surface: smaller y
    loop_upper = gmsh.model.geo.addCurveLoop(
        [l_bottom, l_right1, -l_i3, -l_i2, -l_i1, l_left2]
    )
    s_upper = gmsh.model.geo.addPlaneSurface([loop_upper])

    # Lower surface: larger y
    loop_lower = gmsh.model.geo.addCurveLoop(
        [l_i1, l_i2, l_i3, l_right2, l_top, l_left1]
    )
    s_lower = gmsh.model.geo.addPlaneSurface([loop_lower])

    gmsh.model.geo.synchronize()

    # Physical groups
    pg_upper = gmsh.model.addPhysicalGroup(2, [s_upper])
    gmsh.model.setPhysicalName(2, pg_upper, "upper_layer")

    pg_lower = gmsh.model.addPhysicalGroup(2, [s_lower])
    gmsh.model.setPhysicalName(2, pg_lower, "lower_layer")

    pg_interface = gmsh.model.addPhysicalGroup(1, [l_i1, l_i2, l_i3])
    gmsh.model.setPhysicalName(1, pg_interface, "interface")

    pg_diffr = gmsh.model.addPhysicalGroup(0, [p_k1, p_k2])
    gmsh.model.setPhysicalName(0, pg_diffr, "diffraction_points")

    pg_source = gmsh.model.addPhysicalGroup(0, [p_src])
    gmsh.model.setPhysicalName(0, pg_source, "source")

    # Force source as a mesh node
    gmsh.model.mesh.embed(0, [p_src], 2, s_upper)
    gmsh.model.mesh.embed(0, [p_src], 2, s_lower)

    gmsh.option.setNumber("Mesh.RecombineAll", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    # Nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    xyz = np.array(node_coords, dtype=float).reshape(-1, 3)

    # Triangles
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2)
    tri_type = 2  # 3-node triangle

    tri_elem_tags = None
    tri_nodes = None
    for et, tags_block, nodes_block in zip(elem_types, elem_tags, elem_node_tags):
        if et == tri_type:
            tri_elem_tags = np.array(tags_block, dtype=np.int64)
            tri_nodes = np.array(nodes_block, dtype=np.int64).reshape(-1, 3)
            break

    if tri_nodes is None:
        gmsh.finalize()
        raise RuntimeError("No 3-node triangles found.")

    tag_to_idx = {int(tag): i for i, tag in enumerate(node_tags)}
    tri_conn0 = np.vectorize(lambda t: tag_to_idx[int(t)])(tri_nodes).astype(np.int64)

    tri_vel = np.zeros(len(tri_elem_tags), dtype=float)
    elem_tag_to_row = {int(tag): i for i, tag in enumerate(tri_elem_tags)}

    def assign_velocity(surface_tag, velocity):
        etypes, etags, _ = gmsh.model.mesh.getElements(2, surface_tag)
        for et, tags_block in zip(etypes, etags):
            if et == tri_type:
                for tag in tags_block:
                    tri_vel[elem_tag_to_row[int(tag)]] = float(velocity)

    assign_velocity(s_upper, vup)
    assign_velocity(s_lower, vdown)

    gmsh.finalize()
    return xyz, tri_conn0, tri_vel


# -----------------------------------------------------------------------------
# VTK writer
# -----------------------------------------------------------------------------

def write_vtk_with_point_time(
    filename,
    xyz,
    tri_conn0,
    tri_vel,
    time_at_points,
    title="mesh_velocity_time",
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
        for tri in tri_conn0:
            f.write(f"3 {int(tri[0])} {int(tri[1])} {int(tri[2])}\n")

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


# -----------------------------------------------------------------------------
# Analytical / semi-analytical times
# -----------------------------------------------------------------------------

def interface_y_piecewise(x, xk1, yk1, xk2, yk2):
    """
    Piecewise interface y(x).
    """
    x = np.asarray(x, dtype=float)
    y = np.empty_like(x)

    m = (yk2 - yk1) / (xk2 - xk1)

    mask1 = x <= xk1
    mask2 = (x > xk1) & (x < xk2)
    mask3 = x >= xk2

    y[mask1] = yk1
    y[mask2] = yk1 + m * (x[mask2] - xk1)
    y[mask3] = yk2
    return y


def build_interface_sampling(
    x0, x1, xk1, yk1, xk2, yk2, ds=1.0
):
    """
    Sample I1, I2, I3 at approximately ds = 1 m.
    """
    # I1
    x_I1 = np.arange(x0, xk1 + 0.5 * ds, ds)
    y_I1 = np.full_like(x_I1, yk1, dtype=float)
    I1 = np.column_stack([x_I1, y_I1])

    # I2
    Lr = np.hypot(xk2 - xk1, yk2 - yk1)
    s_I2 = np.arange(0.0, Lr + 0.5 * ds, ds)
    t = s_I2 / Lr
    x_I2 = xk1 + t * (xk2 - xk1)
    y_I2 = yk1 + t * (yk2 - yk1)
    I2 = np.column_stack([x_I2, y_I2])

    # I3
    x_I3 = np.arange(xk2, x1 + 0.5 * ds, ds)
    y_I3 = np.full_like(x_I3, yk2, dtype=float)
    I3 = np.column_stack([x_I3, y_I3])

    return I1, I2, I3, s_I2


def critical_point_on_I1(src, yk1, vup, vdown):
    """
    CA1 on the left horizontal interface for the head-wave from the source.
    """
    sx, sy = src
    ic = np.arcsin(vup / vdown)
    x_ca1 = sx + (yk1 - sy) * np.tan(ic)
    return np.array([x_ca1, yk1], dtype=float), ic


def compute_interface_times(
    src,
    I1,
    I2,
    I3,
    s_I2,
    xk1,
    yk1,
    xk2,
    yk2,
    vup,
    vdown,
):
    """
    Compute the interface times with the separation discussed:

    - direct times on I1 are used for lower-layer/ramp feeding
    - conic on I1 is used only for the upper model
    - ramp times are:
          min( direct-refraction from I1 ,
               T(D1) + propagation along ramp at v2 )
    - I3 times are:
          T(D2) + propagation along I3 at v2

    Returns a dict with all useful arrays/scalars.
    """
    sx, sy = src

    # ---- I1 direct times (for lower layer and ramp feeding)
    d_S_I1 = np.hypot(I1[:, 0] - sx, I1[:, 1] - sy)
    t_I1_direct = d_S_I1 / vup

    # ---- CA1 and conic/head-wave times on I1
    CA1, ic = critical_point_on_I1(src, yk1, vup, vdown)
    t_CA1 = np.hypot(CA1[0] - sx, CA1[1] - sy) / vup
    t_I1_conic = t_CA1 + np.abs(I1[:, 0] - CA1[0]) / vdown
    t_I1_node = np.minimum(t_I1_direct, t_I1_conic)

    # D1 = end of I1
    D1 = np.array([xk1, yk1], dtype=float)
    D2 = np.array([xk2, yk2], dtype=float)

    # ---- Time at D1 used for diffraction branch in the upper model
    # You stated that D1 is linked to the conic wave.
    t_D1 = t_CA1 + np.abs(xk1 - CA1[0]) / vdown

    # ---- Ramp times from direct refraction via I1
    # For each Q on ramp:
    #   min_P [ |S-P|/vup + |P-Q|/vdown ]
    dx = I2[:, 0][:, None] - I1[:, 0][None, :]
    dy = I2[:, 1][:, None] - I1[:, 1][None, :]
    d_I1_I2 = np.hypot(dx, dy)
    t_I2_from_I1 = np.min(t_I1_direct[None, :] + d_I1_I2 / vdown, axis=1)

    # ---- Ramp times from D1 continuation
    t_I2_from_D1 = t_D1 + s_I2 / vdown

    # ---- Final ramp times
    t_I2 = np.minimum(t_I2_from_I1, t_I2_from_D1)

    # ---- D2 time
    t_D2 = t_I2[-1]

    # ---- I3 times
    t_I3 = t_D2 + (I3[:, 0] - xk2) / vdown

    return {
        "CA1": CA1,
        "ic": ic,
        "t_CA1": t_CA1,
        "D1": D1,
        "D2": D2,
        "t_D1": t_D1,
        "t_D2": t_D2,
        "t_I1_direct": t_I1_direct,
        "t_I1_conic": t_I1_conic,
        "t_I1_node": t_I1_node,
        "t_I2": t_I2,
        "t_I3": t_I3,
    }


def upper_branch_I1_conic(nodes_xy, CA1, xk1, yk1, vup, vdown):
    """
    Conic field from I1 in the upper medium.

    We use the standard flat-interface head-wave formula from CA1.
    """
    x = nodes_xy[:, 0]
    y = nodes_xy[:, 1]

    ic = np.arcsin(vup / vdown)
    xoff = np.abs(x - CA1[0])

    crit_off = (yk1 - CA1[1]) * np.tan(ic) + (yk1 - y) * np.tan(ic)
    # CA1 lies on y=yk1, so first term is 0, left here for readability.
    crit_off = (yk1 - y) * np.tan(ic)

    # Flat-interface conic branch
    t = np.full_like(x, np.inf, dtype=float)
    valid = y <= yk1
    valid &= xoff > crit_off

    x2 = (yk1 - y[valid]) * np.tan(ic)
    t[valid] = (yk1 - CA1[1] + yk1 - y[valid]) / (np.cos(ic) * vup)
    # first term above has zero numerator for CA1-on-interface; kept formal
    t[valid] = (yk1 - y[valid]) / (np.cos(ic) * vup)
    t[valid] += (xoff[valid] - x2) / vdown

    # Add travel from source to CA1
    t0 = np.hypot(CA1[0] - 500.0, CA1[1] - 500.0)  # overwritten below by caller style
    return t


def compute_theoretical_times(
    xyz,
    src=(500.0, 500.0),
    x0=0.0,
    width=10000.0,
    xk1=3333.3,
    yk1=1666.7,
    xk2=6666.7,
    yk2=2254.42,
    vup=1000.0,
    vdown=3000.0,
    ds_interface=1.0,
    interface_tol=1.0e-6,
):
    """
    Theoretical-time construction with explicit treatment of interface nodes.

    Node classes:
      - upper nodes: strictly above interface
      - lower nodes: strictly below interface
      - interface nodes: on interface within tolerance

    Interface-node assignment:
      - I1: direct source-to-I1 time
      - I2: ramp time T(Q)
      - I3: T(D2) + propagation along I3 at v2
    """
    x1 = x0 + width
    sx, sy = src

    nodes_xy = xyz[:, :2].copy()

    # -------------------------------------------------------------------------
    # Interface sampling
    # -------------------------------------------------------------------------
    I1, I2, I3, s_I2 = build_interface_sampling(
        x0=x0, x1=x1, xk1=xk1, yk1=yk1, xk2=xk2, yk2=yk2, ds=ds_interface
    )

    iface = compute_interface_times(
        src=src,
        I1=I1,
        I2=I2,
        I3=I3,
        s_I2=s_I2,
        xk1=xk1,
        yk1=yk1,
        xk2=xk2,
        yk2=yk2,
        vup=vup,
        vdown=vdown,
    )

    CA1 = iface["CA1"]
    t_CA1 = iface["t_CA1"]
    D1 = iface["D1"]
    D2 = iface["D2"]
    t_D1 = iface["t_D1"]
    t_D2 = iface["t_D2"]

    t_I1_direct = iface["t_I1_direct"]
    t_I1_node = iface["t_I1_node"]
    t_I2 = iface["t_I2"]
    t_I3 = iface["t_I3"]

    # -------------------------------------------------------------------------
    # Node classification
    # -------------------------------------------------------------------------
    yi = interface_y_piecewise(nodes_xy[:, 0], xk1, yk1, xk2, yk2)
    dyi = nodes_xy[:, 1] - yi

    interface_mask = np.abs(dyi) <= interface_tol
    upper_mask = dyi < -interface_tol
    lower_mask = dyi > interface_tol

    times = np.full(nodes_xy.shape[0], np.inf, dtype=float)

    # -------------------------------------------------------------------------
    # Interface nodes
    # -------------------------------------------------------------------------
    idx_i = np.where(interface_mask)[0]
    if idx_i.size > 0:
        X = nodes_xy[idx_i, 0]
        Y = nodes_xy[idx_i, 1]

        # I1 nodes
        mask_i1 = X <= xk1 + interface_tol

        # I3 nodes
        mask_i3 = X >= xk2 - interface_tol

        # I2 nodes
        mask_i2 = ~(mask_i1 | mask_i3)

        # --- I1: stored nodal first-arrival time on the interface
        if np.any(mask_i1):
            Xi = X[mask_i1]
            times[idx_i[mask_i1]] = np.interp(Xi, I1[:, 0], t_I1_node)

        # --- I2: interpolate ramp time by curvilinear abscissa
        if np.any(mask_i2):
            Xi = X[mask_i2]
            Yi = Y[mask_i2]

            Lr = np.hypot(xk2 - xk1, yk2 - yk1)
            ux = (xk2 - xk1) / Lr
            uy = (yk2 - yk1) / Lr

            s = (Xi - xk1) * ux + (Yi - yk1) * uy
            s = np.clip(s, 0.0, Lr)

            times[idx_i[mask_i2]] = np.interp(s, s_I2, t_I2)

        # --- I3: along-interface time from D2
        if np.any(mask_i3):
            Xi = X[mask_i3]
            times[idx_i[mask_i3]] = t_D2 + (Xi - xk2) / vdown

    # -------------------------------------------------------------------------
    # Upper layer
    # -------------------------------------------------------------------------
    idx_u = np.where(upper_mask)[0]
    if idx_u.size > 0:
        U = nodes_xy[idx_u]

        # Direct wave
        t_dir = np.hypot(U[:, 0] - sx, U[:, 1] - sy) / vup

        # I1 conic branch
        ic = np.arcsin(vup / vdown)
        xoff = np.abs(U[:, 0] - CA1[0])
        crit_off = (yk1 - U[:, 1]) * np.tan(ic)

        t_i1_conic = np.full(U.shape[0], np.inf, dtype=float)
        valid = xoff > crit_off

        if np.any(valid):
            x2 = (yk1 - U[:, 1][valid]) * np.tan(ic)
            t_i1_conic[valid] = t_CA1 + (yk1 - U[:, 1][valid]) / (np.cos(ic) * vup)
            t_i1_conic[valid] += (xoff[valid] - x2) / vdown

        # D1 diffraction
        t_d1 = t_D1 + np.hypot(U[:, 0] - D1[0], U[:, 1] - D1[1]) / vup

        # Ramp re-emission into upper medium
        dx = U[:, 0][:, None] - I2[:, 0][None, :]
        dy = U[:, 1][:, None] - I2[:, 1][None, :]
        d_U_I2 = np.hypot(dx, dy)
        t_i2_up = np.min(t_I2[None, :] + d_U_I2 / vup, axis=1)

        # I3 re-emission into upper medium
        dx = U[:, 0][:, None] - I3[:, 0][None, :]
        dy = U[:, 1][:, None] - I3[:, 1][None, :]
        d_U_I3 = np.hypot(dx, dy)
        t_i3_up = np.min(t_I3[None, :] + d_U_I3 / vup, axis=1)

        # D2 diffraction in the upper medium
        t_d2_up = t_D2 + np.hypot(U[:, 0] - D2[0], U[:, 1] - D2[1]) / vup

        # Final upper time
        # D2 cone naturally competes with I3 / I2 / others
        times[idx_u] = np.minimum.reduce(
            [t_dir, t_i1_conic, t_d1, t_i2_up, t_i3_up, t_d2_up]
        )

    # -------------------------------------------------------------------------
    # Lower layer
    # -------------------------------------------------------------------------
    idx_l = np.where(lower_mask)[0]
    if idx_l.size > 0:
        L = nodes_xy[idx_l]

        # Direct refraction via I1 using ONLY direct time on I1
        dx = L[:, 0][:, None] - I1[:, 0][None, :]
        dy = L[:, 1][:, None] - I1[:, 1][None, :]
        d_L_I1 = np.hypot(dx, dy)
        t_l1 = np.min(t_I1_direct[None, :] + d_L_I1 / vdown, axis=1)

        # Continuation from ramp
        dx = L[:, 0][:, None] - I2[:, 0][None, :]
        dy = L[:, 1][:, None] - I2[:, 1][None, :]
        d_L_I2 = np.hypot(dx, dy)
        t_l2 = np.min(t_I2[None, :] + d_L_I2 / vdown, axis=1)

        # Continuation from I3
        dx = L[:, 0][:, None] - I3[:, 0][None, :]
        dy = L[:, 1][:, None] - I3[:, 1][None, :]
        d_L_I3 = np.hypot(dx, dy)
        t_l3 = np.min(t_I3[None, :] + d_L_I3 / vdown, axis=1)

        # D2 diffraction
        t_d2 = t_D2 + np.hypot(L[:, 0] - D2[0], L[:, 1] - D2[1]) / vdown

        times[idx_l] = np.minimum.reduce([t_l1, t_l2, t_l3, t_d2])

    return times

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    mesh_size = 50.0
    x0 = 0.0
    width = 10000.0
    y0 = 0.0
    height = 5000.0

    xk1 = 3333.3
    yk1 = 1666.7
    xk2 = 6666.7
    yk2 = 2254.42

    src = (500.0, 500.0)

    vup = 1000.0
    vdown = 3000.0

    xyz, tri_conn0, tri_vel = build_complex_two_layer_conformal(
        mesh_size=mesh_size,
        x0=x0,
        width=width,
        y0=y0,
        height=height,
        xk1=xk1,
        yk1=yk1,
        xk2=xk2,
        yk2=yk2,
        src=src,
        vup=vup,
        vdown=vdown,
    )

    theoretical_time = compute_theoretical_times(
        xyz,
        src=src,
        x0=x0,
        width=width,
        xk1=xk1,
        yk1=yk1,
        xk2=xk2,
        yk2=yk2,
        vup=vup,
        vdown=vdown,
        ds_interface=1.0,
    )

    write_vtk_with_point_time(
        "input_velocity_complex.vtk",
        xyz,
        tri_conn0,
        tri_vel,
        theoretical_time,
        title="Complex two-layer mesh with velocity + theoretical point time",
    )

    print("Wrote: input_velocity_complex.vtk")


if __name__ == "__main__":
    main()
