#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 22 13:16:12 2025

@author: sotiriskakaletsis
"""

"""
The SEMFE Heat Transfer Solver
Computational Mechanics

Pre-processing Script
"""

import numpy as np
import xml.etree.ElementTree as ET
from meshpy.triangle import MeshInfo, build
from meshpy.triangle import MeshInfo, build

# ------------------------------------------------------------
# MESH GENERATION (rectangle with concentric rectangular hole)
# ------------------------------------------------------------
def generate_rect_with_hole_mesh(W, H, w, h, max_area=None):
    """
    Generate 2D triangular mesh for:
        Outer rectangle: W × H
        Inner hole:      w × h  (centered)
    """

    # Outer rectangle bounds
    x0, x1 = 0.0, W
    y0, y1 = 0.0, H

    # Hole bounds
    hx0, hx1 = 0.2, 0.6
    hy0, hy1 = 0.2, 0.4

    # Points: outer 4 + inner 4
    points = [
        (x0, y0), (x1, y0), (x1, y1), (x0, y1),     # outer
        (hx0, hy0), (hx1, hy0), (hx1, hy1), (hx0, hy1)  # hole
    ]

    outer_facets = [(0, 1), (1, 2), (2, 3), (3, 0)]
    inner_facets = [(4, 5), (5, 6), (6, 7), (7, 4)]

    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(outer_facets + inner_facets)
    mesh_info.set_holes([(0.4, 0.3)])  # any point inside the hole

    mesh = build(mesh_info, max_volume=max_area)

    nodes = [(i + 1, p[0], p[1]) for i, p in enumerate(mesh.points)]
    elements = [(i + 1, tri[0]+1, tri[1]+1, tri[2]+1)
                for i, tri in enumerate(mesh.elements)]

    return nodes, elements


# ------------------------------------------------------------
# FIND BOUNDARY EDGES + INNER-HOLE NODES
# ------------------------------------------------------------
def find_boundaries(nodes, elements, W, H, w, h):
    node_dict = {nid: (x, y) for nid, x, y in nodes}
    TOL = 1e-8

    x_left, x_right = 0.0, W
    y_bottom, y_top = 0.0, H

    hx0, hx1 = 0.2, 0.6
    hy0, hy1 = 0.2, 0.4

    boundaries = {
        "outer_left": [],
        "outer_right": [],
        "outer_bottom": [],
        "outer_top": [],
        "inner_hole_edges": [],
        "inner_hole_nodes": set(),
        "outer_top_nodes": set(),
    }

    for elem_id, n1, n2, n3 in elements:

        edge_map = {
            1: (n1, n2),
            2: (n2, n3),
            3: (n3, n1),
        }

        for e_id, (a, b) in edge_map.items():
            xa, ya = node_dict[a]
            xb, yb = node_dict[b]

            # --- Outer boundaries ---
            if abs(xa - x_left) < TOL and abs(xb - x_left) < TOL:
                boundaries["outer_left"].append((elem_id, e_id))
            if abs(xa - x_right) < TOL and abs(xb - x_right) < TOL:
                boundaries["outer_right"].append((elem_id, e_id))
            if abs(ya - y_bottom) < TOL and abs(yb - y_bottom) < TOL:
                boundaries["outer_bottom"].append((elem_id, e_id))
            if abs(ya - y_top) < TOL and abs(yb - y_top) < TOL:
                boundaries["outer_top"].append((elem_id, e_id))
                boundaries["outer_top_nodes"].add(a)
                boundaries["outer_top_nodes"].add(b)

            # --- Inner hole edges ---
            if ( (abs(xa-hx1)<TOL and abs(ya)<=hy1) or
                 (abs(xa-hx0)<TOL and abs(ya)<=hy1) or
                 (abs(ya-hy1)<TOL and abs(xa)<=hx1) or
                 (abs(ya-hy0)<TOL and abs(xa)<=hx1)):
                boundaries["inner_hole_nodes"].add(a)
                
            if ( (abs(xb-hx1)<TOL and abs(yb)<=hy1) or
                 (abs(xb-hx0)<TOL and abs(yb)<=hy1) or
                 (abs(yb-hy1)<TOL and abs(xb)<=hx1) or
                 (abs(yb-hy0)<TOL and abs(xb)<=hx1)):
                boundaries["inner_hole_nodes"].add(b)

    return boundaries

def pretty_print_xml(elem):
    """Return a pretty-printed XML string for the Element."""
    import xml.dom.minidom as md
    rough = ET.tostring(elem, encoding="utf-8")
    reparsed = md.parseString(rough)
    return reparsed.toprettyxml(indent="  ", encoding="ISO-8859-1")

# ------------------------------------------------------------
# WRITE SEMFE XML FILE
# ------------------------------------------------------------
def write_semfe_xml(filename, nodes, elements, boundaries):
    root = ET.Element("SEMFE_spec")
    ET.SubElement(root, "Module", {"type": "heat conduction"})

    # Materials
    mats = ET.SubElement(root, "Materials")
    mat = ET.SubElement(mats, "Material", {"id": "1", "name": "Steel"})
    ET.SubElement(mat, "conductivity").text = "1.5"

    # Geometry
    geom = ET.SubElement(root, "Geometry")
    node_blk = ET.SubElement(geom, "Nodes")
    for nid, x, y in nodes:
        ET.SubElement(node_blk, "node", {
            "id": str(nid), "x": str(x), "y": str(y), "z": "0.0"
        })

    elem_blk = ET.SubElement(geom, "Elements",
                             {"type": "tri3", "name": "mesh"})
    for eid, n1, n2, n3 in elements:
        el = ET.SubElement(elem_blk, "elem", {"id": str(eid)})
        el.text = f"{n1} {n2} {n3}"

    # Boundary conditions
    bc = ET.SubElement(root, "BoundaryConditions")

    # ---- (1) Apply constant temperature to INNER HOLE NODES ----
    temperature_blk = ET.SubElement(bc, "Boundary")
    for nid in sorted(boundaries["inner_hole_nodes"]):
        ET.SubElement(temperature_blk, "temperature", {
            "node": str(nid),
            "value": "100.0"       # <<< choose your desired temperature
        })

    temperature_blk = ET.SubElement(bc, "Boundary")
    for nid in sorted(boundaries["outer_top_nodes"]):
        ET.SubElement(temperature_blk, "temperature", {
            "node": str(nid),
            "value": "30.0"       # <<< choose your desired temperature
        })
        
        
    # ---- (2) Apply convection BC on outer-right boundary ----
    conv_blk = ET.SubElement(bc, "Convection")
    for elem_id, edge_id in boundaries["outer_right"]:
        ET.SubElement(conv_blk, "conv", {
            "elem": str(elem_id),
            "edge": str(edge_id),
            "h": "50.0",
            "Tinf": "25.0"
        })
        
    # ---- (3) Apply heatflux on outer-right boundary ----
    htflx_blk = ET.SubElement(bc, "Convection")
    for elem_id, edge_id in boundaries["outer_bottom"]:
        ET.SubElement(htflx_blk, "flux", {
            "elem": str(elem_id),
            "edge": str(edge_id),
            "value": "0.0"
        })

    # Simulation Step
    step = ET.SubElement(root, "Step",
                         {"name": "step1", "type": "steady-state"})
    ET.SubElement(step, "HeatSource")

    ET.ElementTree(root).write(
        filename, encoding="ISO-8859-1", xml_declaration=True
    )
    print(f"Written: {filename}")
    
    # Write pretty XML
    xml_bytes = pretty_print_xml(root)
    with open(filename, "wb") as f:
        f.write(xml_bytes)

    print(f"Written (pretty XML): {filename}")


# ------------------------------------------------------------
# MAIN
# ------------------------------------------------------------
if __name__ == "__main__":
    # Outer rectangle and hole dimensions
    W, H = 0.8, 0.6
    w, h = 0.4, 0.2

    nodesGEOM, elementsGEOM = generate_rect_with_hole_mesh(
            W, H, w, h, max_area=0.00025
        )

    boundariesGEOM = find_boundaries(nodesGEOM, elementsGEOM, W, H, w, h)
    write_semfe_xml("chimney.semfe", nodesGEOM, elementsGEOM, boundariesGEOM)
