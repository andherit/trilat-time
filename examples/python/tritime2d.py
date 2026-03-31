# for calling tritime2d
import ctypes,os
import numpy as np 
import gmsh



_lib_path = os.path.join(os.path.dirname(__file__), '_tritime2d.so')
_lib = ctypes.CDLL(_lib_path)

def tritime2d(px, py, pz, cells, velocity, initial_time, diff_nodes=None, fast=False):
    """
    tritime2d : float64 array shape (Nnodes,), pre-initialised by caller.
                 Set np.inf everywhere except nodes reachable from source.
                 For a node source: traveltime[k] = 0.
                 For an in-cell source: set the cell's nodes to the
                 analytical velocity from source point to each node.
    """
    Nnodes = len(px)
    Ncells = len(velocity)

    px         = np.ascontiguousarray(px,         dtype=np.float64)
    py         = np.ascontiguousarray(py,         dtype=np.float64)
    pz         = np.ascontiguousarray(pz,         dtype=np.float64)
    velocity   = np.ascontiguousarray(velocity,   dtype=np.float64)
    traveltime = np.array(initial_time, dtype=np.float64)  # copy if needed
    cells = np.ascontiguousarray(cells, dtype=np.int32)
    if cells.min() == 0:
        cells = cells + 1   # convert 0-based to 1-based

    if diff_nodes is None or not np.any(diff_nodes):
        diff_nodes = np.zeros(1, dtype=np.int32)
        NdiffNodes = 0
    else:
        # diff_nodes is a mask (length Nnodes): extract 1-based indices of marked nodes
        idx = np.where(np.asarray(diff_nodes) > 0)[0].astype(np.int32) + 1
        diff_nodes = np.ascontiguousarray(idx)
        NdiffNodes = len(diff_nodes)

    c_int_p = ctypes.POINTER(ctypes.c_int)
    c_dbl_p = ctypes.POINTER(ctypes.c_double)

    _lib.tritime2d(
        ctypes.c_int(Nnodes), ctypes.c_int(Ncells),
        px.ctypes.data_as(c_dbl_p),
        py.ctypes.data_as(c_dbl_p),
        pz.ctypes.data_as(c_dbl_p),
        cells.ctypes.data_as(c_int_p),
        velocity.ctypes.data_as(c_dbl_p),
        ctypes.c_int(NdiffNodes),
        diff_nodes.ctypes.data_as(c_int_p),
        ctypes.c_int(int(fast)),
        traveltime.ctypes.data_as(c_dbl_p),
    )
    return traveltime



def extract_mesh():
    """
    Extract nodes and cells from the current gmsh model.
    Returns x, y, z (Nnodes,), cells (Ncells x 3, 0-based), cell_type (Ncells,).
    cell_type holds the gmsh surface tag for each cell — needed for layer-based velocity.
    """
    nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
    idx    = nodeTags.astype(int) - 1
    coords = nodeCoords.reshape(-1, 3)
    nnodes = len(nodeTags)
    x = np.empty(nnodes, dtype=np.float64)
    y = np.empty(nnodes, dtype=np.float64)
    z = np.empty(nnodes, dtype=np.float64)
    x[idx] = coords[:, 0]
    y[idx] = coords[:, 1]
    z[idx] = coords[:, 2]

    cell_list = []
    type_list = []
    for _, s in gmsh.model.getEntities(2):
        _, _, elemNodeTags = gmsh.model.mesh.getElements(2, s)
        n = len(elemNodeTags[0]) // 3
        cell_list.append(elemNodeTags[0].reshape(-1, 3))
        type_list.append(np.full(n, s, dtype=np.int32))
    cells     = np.vstack(cell_list).astype(np.int32) - 1
    cell_type = np.concatenate(type_list)

    return x, y, z, cells, cell_type


def velocity_gradient(cells, y, v_base, v_grad):
    """
    Velocity varies linearly with node depth.
    v_base : initial velocity at surface
    v_grad : gradient of velocity
    """
    return v_base + v_grad * y[cells].mean(axis=1)


def velocity_by_layer(cell_type, layer_v):
    """Velocity set per gmsh surface (layer_v[i] is velocity of surface tag i+1)."""
    return np.asarray(layer_v)[cell_type - 1]


def velocity_constant(ncells, v):
    """Uniform velocity."""
    return np.full(ncells, v, dtype=np.float64)


def init_traveltime(sx, sz, x, y, cells, cell_vel, inf=1.e32):
    """
    Find the cell containing source point (sx, sz), compute initial travel
    times to its nodes, and return a traveltime array (inf everywhere else).
    Replaces the manual cell-search loop.
    """
    # vectorised point-in-triangle test over all cells at once
    v0 = np.column_stack((x[cells[:, 0]], y[cells[:, 0]]))  # (Ncells, 2)
    v1 = np.column_stack((x[cells[:, 1]], y[cells[:, 1]]))
    v2 = np.column_stack((x[cells[:, 2]], y[cells[:, 2]]))
    pt = np.array([sx, sz])

    a = v1 - v0               # (Ncells, 2)
    b = v2 - v0
    c = pt - v0               # broadcasts: (Ncells, 2)

    area = a[:, 0]*b[:, 1] - a[:, 1]*b[:, 0]
    u    = (a[:, 0]*c[:, 1] - a[:, 1]*c[:, 0]) / area
    v    = (c[:, 0]*b[:, 1] - c[:, 1]*b[:, 0]) / area

    hit = np.where((u >= 0) & (v >= 0) & (u + v <= 1))[0]
    if len(hit) == 0:
        raise ValueError(f'Source point ({sx}, {sz}) not found in any cell')
    i = hit[0]

    triangle = np.column_stack((x[cells[i]], y[cells[i]]))  # (3, 2)
    t = np.linalg.norm(pt - triangle, axis=1) / cell_vel[i] # vectorised calc_time

    print(f'Source location [{sx}, {sz}], source is in cell {i}')
    print(f'Initial travel time on starting cell: {t}')

    initial_time = inf * np.ones(len(x), dtype=np.float64)
    near = np.where(t < 1e-5)[0]
    if len(near) > 0:
        initial_time[cells[i, near]] = 0.0
        print('Source located at a single node')
    else:
        initial_time[cells[i]] = t

    return initial_time


def write_vtk(x, y, z, cells, velocity, traveltime, filename):
    """
    Write an unstructured triangular mesh with cell velocity and node traveltime
    fields to a VTK legacy ASCII file (POLYDATA format).

    Parameters
    ----------
    x, y, z    : array-like, node coordinates
    cells      : array-like of shape (ncells, 3), triangle node indices
    velocity   : array-like of length ncells, cell-centred velocity
    traveltime : array-like of length nnodes, node traveltime
    filename   : str, output file path
    """
    nnodes = len(x)
    ncells = len(cells)

    with open(filename, "w") as f:
        # --- header ---
        f.write("# vtk DataFile Version 2.0\n")
        f.write("TravelTime\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        # --- nodes ---
        f.write("POINTS %d double\n" % nnodes)
        for i in range(nnodes):
            f.write("%f %f %f\n" % (x[i], y[i], z[i]))
        f.write("\n")

        # --- connectivity ---
        f.write("POLYGONS %d %d\n" % (ncells, ncells * 4))
        for i in cells:
            f.write("3 %d %d %d\n" % (i[0], i[1], i[2]))
        f.write("\n")

        # --- cell data: velocity ---
        f.write("CELL_DATA %d\n" % ncells)
        f.write("SCALARS velocity double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(ncells):
            f.write("%f\n" % velocity[i])
        f.write("\n")

        # --- point data: traveltime ---
        f.write("POINT_DATA %d\n" % nnodes)
        f.write("SCALARS time double 1\n")
        f.write("LOOKUP_TABLE default\n")
        for i in range(nnodes):
            f.write("%f\n" % traveltime[i])
        f.write("\n")
