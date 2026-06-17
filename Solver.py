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
    area = 0.5*abs(x1*(y2-y3)-x2*(y3-y1)-x3*(y1-y2))# element area

    y23 = y2 - y3
    y31 = y3 - y1
    y12 = y1 - y2
    
    x32 = x3 - x2
    x13 = x1 - x3
    x21 = x2 - x1
    
    # Shape function derivatives (constant over element)
    B = 1/(2*area)*np.array ([[y23, y31, y12], [x32, x13, x21]])
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
                # Assemble global Matrix
                rows.append(i_global)
                cols.append(j_global)
                data.append(Ke[i_local,j_local])
                
    K = sp.coo_matrix((data,(rows,cols)), shape=(nnodes, nnodes)).tocsr()

    return K


def apply_dirichlet(K, f, bc_nodes, bc_values):
    """
    Apply Dirichlet boundary conditions to the global matrix
    bc_nodes: array of node indices
    bc_values: array of prescribed values
    """
    K = K.tolil(copy=True)
    f = f.copy()
    
    for node, val in zip(bc_nodes, bc_values):
        
        #modify K and f accordingly   
        col = K[:, node].toarray().ravel()
        f -= col * val
        
        K[node, :] = 0
        K[:, node] = 0
        K[node, node] = 1
        f[node] = val
        
    return K.tocsr(), f

def apply_dirichlet_penalty(K, f, bc_nodes, bc_values, alpha_factor=1e4):

    K_mod = K.tolil(copy=True)
    f_mod = f.copy()

    Kmax = abs(K_mod).max()
    alpha = Kmax * alpha_factor

    for node, value in zip(bc_nodes, bc_values):
        K_mod[node, node] += alpha
        f_mod[node] += alpha * value

    return K_mod.tocsr(), f_mod


def apply_heat_flux(f, nodes, elems, heat_flux_bcs):
    """
    Apply Neumann (heat flux) BCs to load vector.
    Each BC: (elem_id, edge_id, q)
    """
    f = f.copy()
    for elem_id, edge_id, q in heat_flux_bcs:
        
        print(f"HeatFlux -> element {elem_id+1}, edge {edge_id}, q={q}")
        
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
        fe = q * L * np.array([0.5, 0.5]) # linear edge shape functions
        f[conn[edge_nodes]] += fe
    return f


def apply_convection(K, f, nodes, elems, conv_bcs):
    """
    Apply Robin (convection) BCs to load vector & matrix K.
    Each BC: (elem_id, edge_id, h, Tinf)
    """
    K = K.tolil()
    f = f.copy()
    for elem_id, edge_id, h, Tinf in conv_bcs:
        
        print(f"Convection -> element {elem_id+1}, edge {edge_id}, h={h}, Tinf={Tinf}")
        
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
        #Modify K and F accordingly
      
        Ke = h * L / 6.0 * np.array([[2, 1],
                                     [1, 2]])
        fe = h * L * Tinf / 2.0 * np.array([1.0, 1.0])

        edge_conn = conn[edge_nodes]
        for i_local, i_global in enumerate(edge_conn):
            f[i_global] += fe[i_local]
            for j_local, j_global in enumerate(edge_conn):
                K[i_global, j_global] += Ke[i_local, j_local]

    return K.tocsr(), f

def solve_system(K, f):
    """Solve the linear system Ku=f"""
    u = spla.spsolve(K, f)
    return u