import numpy as np


def build_interface_sampling_hr(x0, x1, xk1, yk1, xk2, yk2, ds=1.0):
    x_I1 = np.arange(x0, xk1 + 0.5 * ds, ds)
    y_I1 = np.full_like(x_I1, yk1, dtype=float)
    I1 = np.column_stack([x_I1, y_I1])

    Lr = np.hypot(xk2 - xk1, yk2 - yk1)
    s_I2 = np.arange(0.0, Lr + 0.5 * ds, ds)
    t = s_I2 / Lr
    x_I2 = xk1 + t * (xk2 - xk1)
    y_I2 = yk1 + t * (yk2 - yk1)
    I2 = np.column_stack([x_I2, y_I2])

    x_I3 = np.arange(xk2, x1 + 0.5 * ds, ds)
    y_I3 = np.full_like(x_I3, yk2, dtype=float)
    I3 = np.column_stack([x_I3, y_I3])

    return I1, I2, I3, s_I2


def critical_point_on_I1(src, yk1, v1, v2):
    sx, sy = src
    theta_c = np.arcsin(v1 / v2)
    x_ca1 = sx + (yk1 - sy) * np.tan(theta_c)
    return np.array([x_ca1, yk1], dtype=float), theta_c


def interp_interface_time_along_s(sample_s, sample_t, node_s):
    return np.interp(node_s, sample_s, sample_t)


def compute_traveltime_bruteforce(
    xyz,
    field_id,
    src=(500.0, 500.0),
    x0=0.0,
    width=10000.0,
    xk1=3333.3,
    yk1=1666.7,
    xk2=6666.7,
    yk2=2254.42,
    v1=1000.0,
    v2=3000.0,
    ds_hr=1.0,
    interface_tol=1.0e-6,
):
    """
    Brute-force traveltime construction using the validated field map.

    Field usage expected:
      0 direct
      1 conic from I1
      2 diffraction from D1
      3 pseudo-conic from I2
      4 transition: min(pseudo-conic from I2, conic from I3)
      5 conic from I3
      6 lower refracted from I1*
      7 shadow from D2

    No field is modified in this stage.
    """
    X = xyz[:, 0]
    Y = xyz[:, 1]
    n = xyz.shape[0]
    sx, sy = src
    x1 = x0 + width

    D1 = np.array([xk1, yk1], dtype=float)
    D2 = np.array([xk2, yk2], dtype=float)
    CA1, theta_c = critical_point_on_I1(src, yk1, v1, v2)
    cos_tc = np.cos(theta_c)
    tan_tc = np.tan(theta_c)

    # Dense interface samplings
    I1_hr, I2_hr, I3_hr, s_I2 = build_interface_sampling_hr(x0, x1, xk1, yk1, xk2, yk2, ds=ds_hr)

    # I1* = from xS to xCA1 on I1
    xmin_i1s = min(sx, CA1[0])
    xmax_i1s = max(sx, CA1[0])
    mask_I1s = (I1_hr[:, 0] >= xmin_i1s - 1e-12) & (I1_hr[:, 0] <= xmax_i1s + 1e-12)
    I1s_hr = I1_hr[mask_I1s]

    # Store source -> Pi1hr time once, reused by f3 and f6
    t_I1s = np.hypot(I1s_hr[:, 0] - sx, I1s_hr[:, 1] - sy) / v1

    # Direct time everywhere in upper fields as a base
    t_direct = np.hypot(X - sx, Y - sy) / v1
    T = np.full(n, np.inf, dtype=float)
    upper_fields = (field_id >= 0) & (field_id <= 5)
    T[upper_fields] = t_direct[upper_fields]

    # Common anchor times
    t_CA1 = np.hypot(CA1[0] - sx, CA1[1] - sy) / v1
    t_D1 = t_CA1 + abs(D1[0] - CA1[0]) / v2

    # ------------------------------------------------------------------
    # f1: conic from I1, compared with direct
    # ------------------------------------------------------------------
    idx = np.where(field_id == 1)[0]
    if idx.size:
        xu = X[idx]
        yu = Y[idx]
        x_ca2 = xu - (yk1 - yu) * tan_tc
        ca1_ca2 = np.abs(x_ca2 - CA1[0])
        ca2_r = (yk1 - yu) / cos_tc
        t_conic_i1 = t_CA1 + ca1_ca2 / v2 + ca2_r / v1
        T[idx] = np.minimum(T[idx], t_conic_i1)

    # ------------------------------------------------------------------
    # f2: diffraction from D1, compared with direct
    # ------------------------------------------------------------------
    idx = np.where(field_id == 2)[0]
    if idx.size:
        t_d1 = t_D1 + np.hypot(X[idx] - D1[0], Y[idx] - D1[1]) / v1
        T[idx] = np.minimum(T[idx], t_d1)

    # ------------------------------------------------------------------
    # I2 dense feeder from I1*
    # T(I2hr) = min( T(I1*) + Pi1*Pi2 / v2 )
    # ------------------------------------------------------------------
    dx = I2_hr[:, 0][:, None] - I1s_hr[:, 0][None, :]
    dy = I2_hr[:, 1][:, None] - I1s_hr[:, 1][None, :]
    d_I1s_I2 = np.hypot(dx, dy)
    t_I2 = np.min(t_I1s[None, :] + d_I1s_I2 / v2, axis=1)
    t_D2 = float(t_I2[-1])

    # Distinct treatment for mesh nodes on I2: linear interpolation along arc length
    idx_i2 = np.where(np.abs((Y - (yk1 + (yk2 - yk1) * (X - xk1) / (xk2 - xk1))) ) <= interface_tol)[0]
    idx_i2 = idx_i2[(X[idx_i2] > xk1 + interface_tol) & (X[idx_i2] < xk2 - interface_tol)]
    if idx_i2.size:
        s_nodes = np.hypot(X[idx_i2] - xk1, Y[idx_i2] - yk1)
        T[idx_i2] = interp_interface_time_along_s(s_I2, t_I2, s_nodes)

    # ------------------------------------------------------------------
    # f3 and f4: pseudo-conic from I2
    # T = min( T(I2hr) + Pi2R / v1 )
    # ------------------------------------------------------------------
    idx34 = np.where((field_id == 3) | (field_id == 4))[0]
    t_pc_i2 = None
    if idx34.size:
        RX = X[idx34]
        RY = Y[idx34]
        dx = RX[:, None] - I2_hr[:, 0][None, :]
        dy = RY[:, None] - I2_hr[:, 1][None, :]
        d_R_I2 = np.hypot(dx, dy)
        t_pc_i2 = np.min(t_I2[None, :] + d_R_I2 / v1, axis=1)

        idx3_mask = field_id[idx34] == 3
        if np.any(idx3_mask):
            T[idx34[idx3_mask]] = t_pc_i2[idx3_mask]

    # ------------------------------------------------------------------
    # f5 and competitor in f4: conic from I3 through D2
    # T = TD2 + D2CA3/v2 + CA3R/v1
    # ------------------------------------------------------------------
    idx45 = np.where((field_id == 4) | (field_id == 5))[0]
    if idx45.size:
        RX = X[idx45]
        RY = Y[idx45]
        x_ca3 = RX - (yk2 - RY) * tan_tc
        d2_ca3 = np.abs(x_ca3 - xk2)
        ca3_r = (yk2 - RY) / cos_tc
        t_conic_i3 = t_D2 + d2_ca3 / v2 + ca3_r / v1

        idx4_mask = field_id[idx45] == 4
        if np.any(idx4_mask):
            # map back to the shared idx34 order if available
            idx4_global = idx45[idx4_mask]
            if t_pc_i2 is None:
                raise RuntimeError("Field 4 present but pseudo-conic from I2 was not computed.")
            t_pc_map = np.empty_like(t_conic_i3[idx4_mask])
            pos34 = {int(g): i for i, g in enumerate(idx34)}
            for k, g in enumerate(idx4_global):
                t_pc_map[k] = t_pc_i2[pos34[int(g)]]
            T[idx4_global] = np.minimum(t_pc_map, t_conic_i3[idx4_mask])

        idx5_mask = field_id[idx45] == 5
        if np.any(idx5_mask):
            T[idx45[idx5_mask]] = t_conic_i3[idx5_mask]

    # ------------------------------------------------------------------
    # f6: lower refracted from I1*
    # T = min( T(I1*) + Pi1*R / v2 )
    # ------------------------------------------------------------------
    idx = np.where(field_id == 6)[0]
    if idx.size:
        RX = X[idx]
        RY = Y[idx]
        dx = RX[:, None] - I1s_hr[:, 0][None, :]
        dy = RY[:, None] - I1s_hr[:, 1][None, :]
        d_R_I1s = np.hypot(dx, dy)
        T[idx] = np.min(t_I1s[None, :] + d_R_I1s / v2, axis=1)

    # ------------------------------------------------------------------
    # f7: shadow cone from D2 in lower medium
    # T = TD2 + D2R / v2
    # ------------------------------------------------------------------
    idx = np.where(field_id == 7)[0]
    if idx.size:
        T[idx] = t_D2 + np.hypot(X[idx] - D2[0], Y[idx] - D2[1]) / v2

    aux = {
        "CA1": CA1,
        "theta_c": theta_c,
        "I1s_hr": I1s_hr,
        "t_I1s": t_I1s,
        "I2_hr": I2_hr,
        "s_I2": s_I2,
        "t_I2": t_I2,
        "t_D1": t_D1,
        "t_D2": t_D2,
    }
    return T, aux
