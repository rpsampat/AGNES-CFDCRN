#!/usr/bin python

import sys
import re
from copy import copy
import pickle
import math

print "Converting from ANSYS Fluent format graph"
# This module has been adapted from Fenicstools; https://github.com/mikaem/fenicstools
# Use regular expressions to identify sections and tokens found in a fluent file
re_dimline = re.compile(r"\(2\s(\d)\)")
re_comment = re.compile(r"\(0\s.*")
re_zone0 = re.compile(r"\(10\s\(0\s(\w+)\s(\w+)\s(\d+)\s(\d+)\)\)")
re_zone = re.compile(r"\(10\s\((\w+)\s(\w+)\s(\w+)\s(\d+)\s(\d)\)(\(|)")
re_face0 = re.compile(r"\(13(\s*)\(0\s+(\w+)\s+(\w+)\s+(0|0 0)\)\)")
re_face = re.compile(r"\(13(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_interface = re.compile(r"\(61(\s*)\((\w+)\s+(\w+)")
re_periodic = re.compile(r"\(18(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_pfaces = re.compile(r"((^\s)|)(\w+)(\s*)(\w+)")
re_cells0 = re.compile(r"\(12(\s*)\(0(\s+)(\w+)(\s+)(\w+)(\s+)(0|0 0)\)\)")
re_cells = re.compile(r"\(12.*\((\w+)\s+(\w+)\s+(\w+)\s+(\d+)\s+(\d+)\)\)")
re_cells2 = re.compile(r"\(12(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\)(\s*)(\(|)")
re_zones = re.compile(r"\((45|39)\s+\((\d+)\s+(\S+)\s+(\S+).*\)\((.*|[0-9]+[\.]*[0-9]*)\)\)")
re_parant = re.compile(r"(\s*\)(\s*)|\s*\)\)(\s*)|\s*\(\s*)")
re_cell_part_id = re.compile(r"\(40(\s*)\((\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s*\)")
# The fluent mesh (the .msh file) is basically stored as a list of nodes, and then a
# list of faces for each zone of the mesh, the interior and the boundaries.

# Declare som maps that will be built when reading in the lists of nodes and faces:
cell_map = {}               # Maps cell id with nodes
cell_face_map = {}          # Maps cell id with faces
face_cell_map = {}          # Maps face id with two cells
face_list = []              # List of faces [[id, 2-4 nodes, 2 connecting cells and type]]
face_map = {}               # For each cell a dictionary with key=face and val=local face number
neighbour = {}              # neighbouring cells
nodes = {}                  # Maps node number to coordinates
pfaces = {}
face_area_list = {}
face_area_except = []
face_node = {}              # maps face id to nodes
node_map = {}               # maps nodes to cells
cell_ax_surf = {}           # maps axial surface area to radial coordinate

# Information about connectivity and boundaries
boundary_cells2face = {}    # {cellnum:{zone_id:[face1,face2..],..},..}

# Some global values
num_cells = {}              # Total number of cells in different zones
zone = {}                   # {id:[cell(12),first index, last index], id:[face(13),fid,lid,type]..}
zones = {}                  # zone information
zone_number_of_faces = {}   # number of faces for each zone
cellzero = []
common = {}


def read_periodic(Nmin, Nmax, ifile):
    """Scan past periodic section"""

    line = ifile.readline()
    readline = False
    if re.search(re_parant, line):  # check for initial paranthesis
        readline = True

    ls = []
    face_index = Nmin
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        string = ''
        lineorig = line
        if re.search(re_parant, line):
            for x in line:
                if x != ')':
                    string = string+x
                # if x==')':
                    # print "last line",i
            line = string
        ln = line.split()
        pfaces[int(ln[0], 16)] = int(ln[1], 16)
        # print "%i:%i"%(int(ln[0],16),int(ln[1],16))


def read_zone_nodes(dim, Nmin, Nmax, ifile):
    """Scan lines for nodes and return in an array."""

    line = ifile.readline()
    readline = False

    if re.search(re_parant, line):  # check for initial paranthesis
        readline = True
        # dummy = lines.pop(0)
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        string = ''
        if re.search(re_parant, line):
            for x in line:
                if x != ')':
                    string = string+x
            line = string
        readline = True
        nodes[i] = [float(x) for x in line.split()]


def face_area():
    """
    calculate face area by splitting polygon into triangles. Assumption that
    nodes in cyclic order in list. Starting node kept as fixed vertex. Other two vertices 
    selected by traversing thorugh list.
    :return: 
    """

    xcoord = []
    ycoord = []
    zcoord = []

    for i in face_node:
        face_area_list[i] = [0, 0, face_cell_map[i]]
            
        """nd=face_node[i]
        if nd==[]:
            face_area_list[i] = [0, 0, face_cell_map[i]]
            continue
        pt1=array(nodes[nd[0]])
        area=0
        for j in range(1,len(nd)-3+2):
            pt2=array(nodes[nd[j]])
            pt3=array(nodes[nd[j+1]])
            a=linalg.norm(pt1-pt2)
            b=linalg.norm(pt2-pt3)
            c=linalg.norm(pt1-pt3)
            s=(a+b+c)/2
            area+=math.sqrt(s*(s-a)*(s-b)*(s-c))
        v1 = subtract(pt1, pt3)
        v2 = subtract(pt3, pt2)
        v = cross(v1, v2)
        face_area_list[i] = [area, v, face_cell_map[i]]"""

    for j in face_area_except:
        try:
            del face_area_list[j]
        except KeyError:
            continue
                

def read_faces(zone_id, Nmin, Nmax, bc_type, face, ifile):
    """Read all faces and create cell_face_map + some boundary maps."""

    line = ifile.readline()
    readline = False
    if re.search(re_parant, line):  # check for initial paranthesis
        readline = True

    ls = []
    face_index = Nmin
    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        string = ''
        lineorig = line
        if re.search(re_parant, line):
            for x in line:
                if x != ')':
                    string = string+x
                # if x==')':
                    # print "last line",i
            line = string
        ln = line.split()
        if face == 0:
            nd = int(ln[0], 16)  # Number of nodes
            nds = [int(x, 16) for x in ln[1:(nd + 1)]]
            cells = [int(x, 16) for x in ln[(nd + 1):]]
            face_index = i
        else:
            nd = face
            nds = [int(x, 16) for x in ln[:nd]]
            cells = [int(x, 16) for x in ln[nd:]]

        face_node[i] = nds
        face_cell_map[i] = cells
        
        for c in cells:
            if c > 0:
                try:
                    cell_map[c] = list(set(cell_map[c] + nds))
                    neighbour[c]['adjacent'] = list(set(neighbour[c]['adjacent'] + cells))
                    neighbour[c]['faces'].append(i)
                except:
                    cell_map[c] = copy(nds)
                    neighbour[c] = {}
                    neighbour[c]['adjacent'] = copy(cells)
                    neighbour[c]['faces'] = [i]
                    
                for n in nds:
                    try:
                        node_map[n].append(c)
                    except KeyError:
                        node_map[n] = [c]
                    
                # neighbour[c]['adjacent'].remove(c)
                neighbour[c]['color'] = 'white'
                neighbour[c]['distance'] = 'inf'
                neighbour[c]['predecessor'] = ''
                neighbour[c]['cells'] = [c]

                """if c in cellzero:
                    try:
                        common[zone_id].append(c)
                    except KeyError:
                        common[zone_id]=[c]
            if c > 0 and nds==[]:
                cellzero.append(c)"""


def interface(Nmin, Nmax, ifile):
    line = ifile.readline()
    readline = False
    if re.search(re_parant, line):  # check for initial paranthesis
        readline = True

    for i in range(Nmin, Nmax + 1):
        if readline:
            line = ifile.readline()
        readline = True
        string = ''
        if re.search(re_parant, line):
            for x in line:
                if x != ')':
                    string = string+x
                # if x==')':
                    # print "last line",i
            line = string
        ln = line.split()
        
        try:
            f1 = int(ln[0], 16)
            f2 = int(ln[1], 16)
        except:
            continue
        if f1 == 0 or f2 == 0:
            continue
        c1 = face_cell_map[f1]
        c2 = face_cell_map[f2]
        c = face_cell_map[i]
        # neighbour[c[0]]['adjacent'].remove(c[1])
        # neighbour[c[1]]['adjacent'].remove(c[0])
        
        # neighbour[c1[0]]['adjacent'].append(c2[0])
        # neighbour[c1[0]]['faces'].append(f2)
        
        # neighbour[c2[0]]['adjacent'].append(c1[0])
        # neighbour[c2[0]]['faces'].append(f1)
        # face_area_except.append(f1)
        # face_area_except.append(f2)
        # print face_cell_map[i]
        # print "c1=",c1
        # print "c2=",c2
        # face_cell_map[i] = [c1[0],c2[0]]
        # face_node[i] = []


def scan_fluent_mesh(ifile):
    """Scan fluent mesh and generate numerous maps."""

    # Warning! Not yet tested for multiple interior zones
    dim = 0
    one = 0
    num_faces = 0
    while 1:
        line = ifile.readline()
        if len(line) == 0:
            print 'Finished reading file\n'
            break

        if dim == 0:  # Dimension usually comes first
            a = re.search(re_dimline, line)
            if a:
                print 'Reading dimensions\n'
                dim = int(a.group(1))
                print 'Mesh is ' + str(dim) + 'D\n'
                continue

        if one == 0:  # The total number of nodes
            a = re.search(re_zone0, line)
            if a:
                print 'Reading zone info\n'
                one, num_vertices, dummy1, dummy2 = int(a.group(1)), int(a.group(2), 16), int(a.group(3), 16), \
                                                    int(a.group(4))
                continue

        a = re.search(re_zone, line)  # Nodes
        if a:
            zone_id, first_id, last_id, dummy1, dummy2 = int(a.group(1), 16), int(a.group(2), 16), int(a.group(3), 16),\
                                                         int(a.group(4)), int(a.group(5))
            print 'Reading ', last_id - first_id + 1, ' nodes in zone ', zone_id + 1, '\n'
            read_zone_nodes(dim, first_id, last_id, ifile)
            continue

        a = re.search(re_zones, line)  # Zone info
        if a:
            print 'Reading zone ', line
            dummy, zone_id, zone_type, zone_name, radius = int(a.group(1)), int(a.group(2)),  a.group(3), a.group(4), \
                                                           a.group(5)
            zones[zone_id] = [zone_type, zone_name, radius]
            continue

        a = re.search(re_cells0, line)  # Get total number of cells/elements
        if a:
            print 'Reading cell info ', line
            first_id, tot_num_cells = int(a.group(3), 16), int(a.group(5), 16)
            continue

        a = re.search(re_cells, line)  # Get the cell info.
        if a:
            zone_id, first_id, last_id, bc_type, element_type = int(a.group(1), 16), int(a.group(2), 16), \
                                                                int(a.group(3), 16), int(a.group(4), 16), \
                                                                int(a.group(5), 16)
            print 'Reading ', last_id - first_id + 1, ' cells in zone ', zone_id, '\n'
            if last_id == 0:
                raise TypeError("Zero elements!")
            num_cells[zone_id] = [first_id, last_id, bc_type, element_type]
            zone[zone_id] = [12, first_id, last_id]
            continue

        a = re.search(re_cells2, line)  # Get the cell info.
        if a:
            print line
            raise TypeError("Wrong cell type. Can only handle one single cell type")

        a = re.search(re_face0, line)
        if a:
            print 'Reading total number of faces\n', line
            num_faces = int(a.group(3), 16)
            continue

        a = re.search(re_face, line)
        if a:
            print 'Reading faces ', line
            print line
            zone_id, first_id, last_id, bc_type, face_type = int(a.group(2), 16), int(a.group(3), 16), \
                                                             int(a.group(4), 16), int(a.group(5), 16), \
                                                             int(a.group(6), 16)

            if bc_type == 12:
                print "Periodic zone"
                print zone_id
                
            # print zone_id
            # print first_id
            # print last_id
            read_faces(zone_id, first_id, last_id, bc_type, face_type, ifile)
            zone_number_of_faces[zone_id] = last_id - first_id + 1
            zone[zone_id] = [13, first_id, last_id, bc_type]
            continue

        a = re.search(re_interface, line)
        if a:
            print "Reading Interface", line
            Nmin, Nmax = int(a.group(2), 16), int(a.group(3), 16)
            interface(Nmin, Nmax, ifile)

        a = re.search(re_cell_part_id, line)
        if a:
            print 'Reading cell partition id ', line
            line = ifile.readline()
            if re.search(re_parant, line):  # check for initial paranthesis
                while 1:
                    # read through lines
                    line = ifile.readline()
                    if re.search(re_parant, line):  # check for closing paranthesis
                        break
            continue

        a = re.search(re_periodic, line)
        if a:
            print 'Reading periodic connectivity\n', line
            Nmin, Nmax = int(a.group(2), 16), int(a.group(3), 16)
            read_periodic(Nmin, Nmax, ifile)
            continue

        # print 'Line = ',line
        if any([re.search(st, line) for st in (re_parant, re_comment)]) or not line.strip():
            # print [re.search(st, line) for st in (re_parant, re_comment)]
            continue

        # Should not make it here
        # print 'Line = ',line
        # raise IOError('Something went wrong reading fluent mesh.')


def convert(fluentmesh):
    """Converts a fluent mesh to a mesh format that can be used by FEniCS.

         fluentmesh = fluent mesh (*.msh file)

    """

    ofilename = fluentmesh[:-4]
    ifile = open(fluentmesh, "r")
    scan_fluent_mesh(ifile)
    ifile.close()
    # print "Common:\n",common
    print "Adding cells with common nodes to graph"
    # adding cells with common nodes to graph
    for cell in cell_map:
        neighbour[cell]['nodes'] = cell_map[cell]
        for node in cell_map[cell]:
            neighbour[cell]['adjacent'] += node_map[node]
            # neighbour[cell]['faces']+=neighbour[node_map[node]]['faces']
        neighbour[cell]['adjacent'] = set(neighbour[cell]['adjacent'])
        # neighbour[cell]['faces'] = set(neighbour[cell]['faces'])
        # Removing cell self id from graph
        while 1:
            try:
                neighbour[cell]['adjacent'].remove(cell)                
            except:
                break
        # Removing zeros from adjacent
       
        while 1:
            try:
                neighbour[cell]['adjacent'].remove(0)                
            except:
                break
        # Removing zeros from faces
        neighbour[cell]['faces'] = set(neighbour[cell]['faces'])
        while 1:
            try:
                neighbour[cell]['faces'].remove(0)                
            except:
                break

    print "Writing Mesh"
    with open('mesh.pkl', 'wb') as file:
        pickle.dump(neighbour, file, pickle.HIGHEST_PROTOCOL)
    print "Writing Zones"
    with open('zone.pkl', 'wb') as file:
        pickle.dump(zone, file, pickle.HIGHEST_PROTOCOL)
    print "Writing Boundary cells"
    with open('bc2f.pkl', 'wb') as file:
        pickle.dump(boundary_cells2face, file, pickle.HIGHEST_PROTOCOL)
    print "Writing periodic faces"
    with open('pfaces.pkl', 'wb') as file:
        pickle.dump(pfaces, file, pickle.HIGHEST_PROTOCOL)
    print "Creating face area list"
    face_area()
    try:
        for face in face_node.keys():  # adding the surface area of faces to face_area_list TODO: make 3D appropriate
            Delta_x = nodes[face_node[face][0]][0] - nodes[face_node[face][1]][0]
            Delta_y = nodes[face_node[face][0]][1] - nodes[face_node[face][1]][1]
            face_length = math.sqrt(Delta_x**2 + Delta_y**2)
            face_surf = face_length * (min(nodes[face_node[face][0]][1], nodes[face_node[face][1]][1]) + 0.5 * abs(Delta_y))
            face_area_list[face].append(face_surf)
            if Delta_x < Delta_y:
                cell_ax_surf[float(str(0.5*(nodes[face_node[face][0]][1] + nodes[face_node[face][1]][1]))[:6])] = face_surf
    except:
        print ""
    print "Writing face area list"
    with open('facearea.pkl', 'wb') as file:
        pickle.dump(face_area_list, file, pickle.HIGHEST_PROTOCOL)
    try:
        with open('cellaxsurf.pkl', 'wb') as file:
            pickle.dump(cell_ax_surf, file, pickle.HIGHEST_PROTOCOL)
    except:
        print ""
    print "Finished Writing Files"
    return 0


if __name__ == '__main__':
    neighb = convert("CFD_CRN.cas")
