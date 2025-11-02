#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The SEMFE Heat Transfer Solver
Computational Mechanics

Solver Script
"""

# --------------------------
# File: fem.py
# --------------------------
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import numpy as np


def element_stiffness_triangle(node_coords, k=1.0):
    """
    Linear triangular element stiffness for steady-state conduction (Poisson equation)
    node_coords: (3,2) or (3,3) array of node coordinates
    returns 3x3 element stiffness matrix
    """
    x1, y1 = node_coords[0, 0], node_coords[0, 1]
    x2, y2 = node_coords[1, 0], node_coords[1, 1]
    x3, y3 = node_coords[2, 0], node_coords[2, 1]
    area = []# element area

    # Shape function derivatives (constant over element)
    B = [] #matrix B
    Ke = k * area * (B.T @ B)
    return Ke


def assemble_global(nodes, elems, k=1.0):
    """
    Assemble global stiffness matrix for triangular mesh
    nodes: Nx2 or Nx3 array
    elems: Mx3 array of node indices (0-based)
    k: thermal conductivity
    returns: K (sparse CSR matrix)
    """
    nnodes = nodes.shape[0]
    nelems = elems.shape[0]
    rows = []
    cols = []
    data = []

    for e in range(nelems):
        conn = elems[e]
        coords = nodes[conn, :2]  # take x,y only
        Ke = element_stiffness_triangle(coords, k=k)
        for i_local, i_global in enumerate(conn):
            for j_local, j_global in enumerate(conn):
                # Assemble global Matrix K
                K=[]

    K = [] 
    return K


def apply_dirichlet(K, f, bc_nodes, bc_values):
    """
    Apply Dirichlet boundary conditions to the global matrix
    bc_nodes: array of node indices
    bc_values: array of prescribed values
    """
    K = K.tocsr().copy()
    f = f.copy()
    for node, val in zip(np.atleast_1d(bc_nodes), np.atleast_1d(bc_values)):
        #modify K and f accordingly    
        K[node, node] = []
        f[node] = []
    return K, f

def apply_heat_flux(f, nodes, elems, heat_flux_bcs):
    """
    Apply Neumann (heat flux) BCs to load vector.
    Each BC: (elem_id, edge_id, q)
    """
    for elem_id, edge_id, q in heat_flux_bcs:
        conn = elems[elem_id]
        coords = nodes[conn, :2]
        edge_nodes = {
            1: [0,1],
            2: [1,2],
            3: [2,0]
        }[edge_id]
        n1, n2 = edge_nodes
        x1, y1 = coords[n1]
        x2, y2 = coords[n2]
        L = np.hypot(x2-x1, y2-y1)
        fe = [] # linear edge shape functions
        f[conn[edge_nodes]] += fe
    return f


def apply_convection(K, f, nodes, elems, conv_bcs):
    """
    Apply Robin (convection) BCs to load vector & matrix K.
    Each BC: (elem_id, edge_id, h, Tinf)
    """
    for elem_id, edge_id, h, Tinf in conv_bcs:
        conn = elems[elem_id]
        coords = nodes[conn, :2]
        edge_nodes = {
            1: [0,1],
            2: [1,2],
            3: [2,0]
        }[edge_id]
        n1, n2 = edge_nodes
        x1, y1 = coords[n1]
        x2, y2 = coords[n2]
        #Modify K and F accordingly
        K = []
        f = []
    return K, f


def solve_system(K, f):
    """Solve the linear system Ku=f"""
    u = spla.spsolve(K, f)
    return u
