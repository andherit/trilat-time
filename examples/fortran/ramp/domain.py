import gmsh
import numpy as np
from domain_time_bruteforce import compute_traveltime_bruteforce


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
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("complex_two_layer_branch_map")

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
    p_k1 = gmsh.model.geo.addPoint(xk1, yk1, 0.0, mesh_size)  # D1
    p_k2 = gmsh.model.geo.addPoint(xk2, yk2, 0.0, mesh_size)  # D2
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

def write_vtk_with_branch_id(
    filename,
    xyz,
    tri_conn0,
    tri_vel,
    branch_id,
    traveltime,
    title="mesh_velocity_branch_id",
):
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
        for tri in tri_conn0:
            f.write(f"3 {int(tri[0])} {int(tri[1])} {int(tri[2])}\n")

        f.write(f"\nCELL_DATA {ntri}\n")
        f.write("SCALARS velocity float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for v in tri_vel:
            f.write(f"{float(v):12.5f}\n")

        f.write(f"\nPOINT_DATA {npts}\n")
        f.write("SCALARS branch_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for b in branch_id:
            f.write(f"{int(b)}\n")

        f.write("SCALARS traveltime float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for t in traveltime:
            f.write(f"{float(t):12.5f}\n")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

def interface_y_piecewise(x, xk1, yk1, xk2, yk2):
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


def build_interface_sampling(x0, x1, xk1, yk1, xk2, yk2, ds=1.0):
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
    sx, sy = src
    ic = np.arcsin(vup / vdown)
    x_ca1 = sx + (yk1 - sy) * np.tan(ic)
    return np.array([x_ca1, yk1], dtype=float), ic


def compute_ramp_and_aux_times(
    src,
    I1,
    I2,
    I3,
    xk1,
    yk1,
    xk2,
    yk2,
    vup,
    vdown,
):
    sx, sy = src

    CA1, ic = critical_point_on_I1(src, yk1, vup, vdown)
    D1 = np.array([xk1, yk1], dtype=float)
    D2 = np.array([xk2, yk2], dtype=float)

    # I1* = from source vertical offset to CA1
    x_os = sx
    x_ca = CA1[0]

    mask_I1_star = (I1[:, 0] >= min(x_os, x_ca) - 1e-12) & (I1[:, 0] <= max(x_os, x_ca) + 1e-12)
    I1_star = I1[mask_I1_star]

    # Direct source -> I1*
    t_S_I1_star = np.hypot(I1_star[:, 0] - sx, I1_star[:, 1] - sy) / vup

    # Branch 3 on I2: min over I1* of S-P1 + P1-Q
    dx = I2[:, 0][:, None] - I1_star[:, 0][None, :]
    dy = I2[:, 1][:, None] - I1_star[:, 1][None, :]
    d_I1s_I2 = np.hypot(dx, dy)
    t_I2_b3 = np.min(t_S_I1_star[None, :] + d_I1s_I2 / vdown, axis=1)

    # Time at D2 from branch 3
    t_D2 = t_I2_b3[-1]

    # Branch 6 feeder on I3
    t_I3_b6 = t_D2 + (I3[:, 0] - xk2) / vdown

    # Time at D1 from I1 conic
    t_CA1 = np.hypot(CA1[0] - sx, CA1[1] - sy) / vup
    t_D1 = t_CA1 + abs(xk1 - CA1[0]) / vdown

    return {
        "CA1": CA1,
        "ic": ic,
        "D1": D1,
        "D2": D2,
        "I1_star": I1_star,
        "t_S_I1_star": t_S_I1_star,
        "t_I2_b3": t_I2_b3,
        "t_I3_b6": t_I3_b6,
        "t_D1": t_D1,
        "t_D2": t_D2,
    }


def physical_angle_from_xy_vector(vec_xy):
    # Convert to standard math coordinates (y up) before atan2
    return np.arctan2(-vec_xy[..., 1], vec_xy[..., 0])


def compute_branch7_wedge_mask(U, D2, left_ray_xy, right_ray_xy):
    """
    U: (N,2) upper nodes
    D2: (2,)
    left_ray_xy, right_ray_xy: unit vectors in model coordinates (x right, y down)

    Returns mask for the wedge above D2 between the two rays.
    """
    V = U - D2[None, :]
    above = V[:, 1] < 0.0

    a = physical_angle_from_xy_vector(V)
    a_left = physical_angle_from_xy_vector(left_ray_xy)
    a_right = physical_angle_from_xy_vector(right_ray_xy)

    a_min = min(a_left, a_right)
    a_max = max(a_left, a_right)

    return above & (a >= a_min) & (a <= a_max)


def compute_left_boundary_ray_from_D2(src, I1_star, t_S_I1_star, D2, xk1, yk1, xk2, yk2, vup, vdown):
    """
    Left boundary of branch 7:
    transmitted ray at D2 that satisfies Snell law across the ramp
    for the incident lower-layer ray arriving at D2 from the optimal I1* -> D2 path.
    """
    # Find optimal P* on I1* for D2
    d_to_D2 = np.hypot(I1_star[:, 0] - D2[0], I1_star[:, 1] - D2[1]) / vdown
    j = np.argmin(t_S_I1_star + d_to_D2)
    Pstar = I1_star[j]

    # Incident direction in lower layer toward D2
    u_in = D2 - Pstar
    u_in = u_in / np.hypot(u_in[0], u_in[1])

    # Ramp tangent and upward normal
    t = np.array([xk2 - xk1, yk2 - yk1], dtype=float)
    t = t / np.hypot(t[0], t[1])

    n_up = np.array([t[1], -t[0]], dtype=float)   # points toward upper medium
    n_up = n_up / np.hypot(n_up[0], n_up[1])

    # Tangential component preserved in Snell construction
    tau_down = np.dot(u_in, t)
    tau_up = (vup / vdown) * tau_down
    tau_up = np.clip(tau_up, -1.0, 1.0)

    beta_up = np.sqrt(max(0.0, 1.0 - tau_up * tau_up))
    u_tr = tau_up * t + beta_up * n_up
    u_tr = u_tr / np.hypot(u_tr[0], u_tr[1])

    return u_tr


# -----------------------------------------------------------------------------
# Branch computation
# -----------------------------------------------------------------------------

def compute_field_id(
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
    nodes_xy = xyz[:, :2].copy()
    CA1, _ = critical_point_on_I1(src, yk1, vup, vdown)
    D1 = np.array([xk1, yk1], dtype=float)
    D2 = np.array([xk2, yk2], dtype=float)

    theta_c = np.arcsin(vup / vdown)
    theta_r = np.arctan2(yk2 - yk1, xk2 - xk1)

    m_crit = 1.0 / np.tan(theta_c)
    aI2 = (yk2 - yk1) / (xk2 - xk1)
    bI2 = yk1 - aI2 * xk1

    m_CA1 = -m_crit
    m_D1I1 = -m_crit
    m_D1I2 = np.tan(np.arctan(m_D1I1) + theta_r)

    # Upward-right rays from D2 must have negative slope since y is positive downward
    m_D2I3 = -m_crit
    m_D2I2 = np.tan(np.arctan(m_D2I3) + theta_r)

    def line_y(x, xref, yref, slope):
        return slope * (x - xref) + yref

    yi = interface_y_piecewise(nodes_xy[:, 0], xk1, yk1, xk2, yk2)
    dyi = nodes_xy[:, 1] - yi

    upper_mask = dyi < -interface_tol
    lower_mask = dyi > interface_tol
    interface_mask = np.abs(dyi) <= interface_tol

    field_id = np.full(nodes_xy.shape[0], -1, dtype=int)

    # -------------------------------------------------------------------------
    # Interface nodes
    # -------------------------------------------------------------------------
    idx_i = np.where(interface_mask)[0]
    if idx_i.size > 0:
        Xi = nodes_xy[idx_i, 0]

        # Explicit split of the 3 interface segments
        mask_i1 = Xi <= xk1 + interface_tol
        mask_i2 = (Xi > xk1 + interface_tol) & (Xi < xk2 - interface_tol)
        mask_i3 = Xi >= xk2 - interface_tol

        if np.any(mask_i1):
            Xi1 = Xi[mask_i1]
            field_id[idx_i[mask_i1]] = np.where(Xi1 < CA1[0], 0, 1)

        if np.any(mask_i2):
            field_id[idx_i[mask_i2]] = 6

        if np.any(mask_i3):
            field_id[idx_i[mask_i3]] = 7

    # -------------------------------------------------------------------------
    # Lower nodes
    # -------------------------------------------------------------------------
    idx_l = np.where(lower_mask)[0]
    if idx_l.size > 0:
        L = nodes_xy[idx_l]
        y_ramp = aI2 * L[:, 0] + bI2
        field_id[idx_l] = np.where(L[:, 1] > y_ramp, 6, 7)

    # -------------------------------------------------------------------------
    # Upper nodes
    # -------------------------------------------------------------------------
    idx_u = np.where(upper_mask)[0]
    if idx_u.size > 0:
        U = nodes_xy[idx_u]
        X = U[:, 0]
        Y = U[:, 1]

        field_id[idx_u] = 0

        y_CA1 = line_y(X, CA1[0], CA1[1], m_CA1)
        y_D1I1 = line_y(X, D1[0], D1[1], m_D1I1)
        y_D1I2 = line_y(X, D1[0], D1[1], m_D1I2)
        y_D2I3 = line_y(X, D2[0], D2[1], m_D2I3)
        y_D2I2 = line_y(X, D2[0], D2[1], m_D2I2)

        mask_f1 = upper_mask[idx_u] & (Y > y_CA1) & (Y < y_D1I1)
        mask_f2 = upper_mask[idx_u] & (Y > y_D1I1) & (Y < y_D1I2)
        mask_f3 = upper_mask[idx_u] & (Y > y_D1I2) & (Y < y_D2I3)
        y_D2_low = np.minimum(y_D2I3, y_D2I2)
        y_D2_high = np.maximum(y_D2I3, y_D2I2)

        mask_f4 = upper_mask[idx_u] & (Y > y_D2_low) & (Y < y_D2_high)
        mask_f5 = upper_mask[idx_u] & (Y > y_D2_high)

        field_id[idx_u[mask_f1]] = 1
        field_id[idx_u[mask_f2]] = 2
        field_id[idx_u[mask_f3]] = 3
        field_id[idx_u[mask_f4]] = 4
        field_id[idx_u[mask_f5]] = 5

    return field_id


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

    branch_id = compute_field_id(
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
        interface_tol=1.0e-6,
    )

    T, aux = compute_traveltime_bruteforce(
        xyz,
        field_id=branch_id,
        src=src,
        x0=x0,
        width=width,
        xk1=xk1,
        yk1=yk1,
        xk2=xk2,
        yk2=yk2,
        v1=vup,
        v2=vdown,
        ds=1.0,
        interface_tol=1.0e-6,
    )

    write_vtk_with_branch_id(
        "input_velocity_complex_branch_id.vtk",
        xyz,
        tri_conn0,
        tri_vel,
        branch_id,
        T,
        title="Complex two-layer mesh with velocity + branch_id",
    )

    print("Wrote: input_velocity_complex_branch_id.vtk")


if __name__ == "__main__":
    main()
