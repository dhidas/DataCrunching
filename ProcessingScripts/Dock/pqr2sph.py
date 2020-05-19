#!/usr/bin/env python3

from array import *
from math  import *

def parse_arguments():
    '''
    Parse the arguments.
    '''
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("pqr",help="pqr file containing fpocket spheres defining pocket")
    args=parser.parse_args()
    return args

def parse_pqr_line(line):
    '''
    Take a line from a PQR file and return a tuple containing
    (xcoord,ycoord,zcoord,radius) where the x, y, and z coordinates
    are the coordinates of the alpha-sphere (see fpocket), and
    the radius is the radius of the alpha-sphere.
    '''
    length=len(line)
    if length < 71:
      print("line too short: ",line)
      return None
    xcoord=float(line[31:39])
    ycoord=float(line[39:47])
    zcoord=float(line[47:55])
    radius=float(line[67:])
    return (xcoord,ycoord,zcoord,radius)

def write_sph_line(data):
    '''
    Take the (xcoord,ycoord,zcoord,radius) tuple and generate a sphere
    line, and return that line.
    '''
    (xcoord,ycoord,zcoord,radius)=data
    line=f"   63{xcoord:10.5f}{ycoord:10.5f}{zcoord:10.5f}{radius:8.3f}   92 0  0"
    return line

def write_pdb_line(data):
    '''
    Take the (xcoord,ycoord,zcoord,radius) tuple and generate a sphere
    line, and return that line.
    '''
    (xcoord,ycoord,zcoord,radius)=data
    line=f"ATOM      5    C STP A   1    {xcoord:8.3f}{ycoord:8.3f}{zcoord:8.3f}  0.00  0.00           C"
    return line

def distance(pnt_i,pnt_j):
    '''
    Calculate the Euclidian distance between two points.
    '''
    (xi,yi,zi,ri) = pnt_i
    (xj,yj,zj,rj) = pnt_j
    dx = xi-xj
    dy = yi-yj
    dz = zi-zj
    r = sqrt(dx*dx+dy*dy+dz*dz)
    return r

def prune_points(data_in,num_out):
    '''
    Go through the list points in the input, compute their relative
    distances, and remove points that are close together until
    only num_out points remain.
    Return the list of remaining points.
    '''
    data_out = []
    del_list = []
    huge = 1.0e100
    num_in = len(data_in)
    num_del = 0
    R = [ [ 0.0 for _ in range(num_in) ] for _ in range(num_in) ]
    for i in range(num_in):
        data_i = data_in[i]
        for j in range(num_in):
            data_j = data_in[j]
            R[i][j] = distance(data_i,data_j)
        R[i][i] = huge
    num_left = num_in - num_del
    while (num_left > num_out):
        min_d = R[0][0]
        min_i = 0
        min_j = 0
        for i in range(num_in): 
            if i in del_list:
                continue
            for j in range(num_in): 
                if j in del_list:
                    continue
                if (R[i][j] < min_d):
                    min_d = R[i][j]
                    min_i = i
                    min_j = j
        min_di = huge
        min_dj = huge
        min_ki = -1
        min_kj = -1
        for k in range(num_in):
            if (k != min_j and R[min_i][k] < min_di):
                min_di = R[min_i][k]
                min_ki = k
            if (k != min_i and R[min_j][k] < min_dj):
                min_dj = R[min_j][k]
                min_kj = k
        if (min_di <= min_dj):
            del_list.append(min_i)
            for k in range(num_in):
                R[min_i][k] = huge
                R[k][min_i] = huge
        else:
            del_list.append(min_j)
            for k in range(num_in):
                R[min_j][k] = huge
                R[k][min_j] = huge
        num_del = len(del_list)
        num_left = num_in - num_del
    for i in range(num_in):
        if i in del_list:
            continue
        data_out.append(data_in[i])
    return data_out


def process_lines(list_in):
    '''
    Go through all lines and convert each one to the spheres format, add it to the result,
    and return the list of output lines.
    '''
    data_in=[]
    list_out=[]
    length=len(list_in)
    list_out.append(f"DOCK spheres within 10.0 ang of ligands")
    list_out.append(f"cluster     1   number of spheres in cluster  {length:4d}")
    for line_in in list_in:
        keyw=line_in[0:6]
        if keyw == b"ATOM  " or keyw == b"HETATM":
            data=parse_pqr_line(line_in)
            data_in.append(data)
    data_out = prune_points(data_in,100)
    for data in data_out:
        #line_out=write_pdb_line(data)
        line_out=write_sph_line(data)
        list_out.append(line_out)
    return list_out

if __name__ == "__main__":
    args=parse_arguments()
    fobj=open(args.pqr,"rb")
    pqr_in=fobj.readlines()
    fobj.close()
    sph_out=process_lines(pqr_in)
    for line_out in sph_out:
        print(line_out)
