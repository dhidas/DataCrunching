#!/usr/bin/env python3
'''
Compute the grid box for a receptor pocket

In drug-screening Fpocket may be used to find a binding pocket in a protein.
A pocket is given by a collection of alpha-spheres. Each alpha-sphere is 
defined as a sphere that touches 4 atoms, demarcating the opening of 
the pocket. Each sphere is given by its position and radius in a
.pqr file.

The grid program builds a docking grid on which the potentials for
docking are expressed. The grid is defined in a rectangular box that
encompasses the binding pocket. The box is given by the positions of its
corners.

This script takes a collection of alpha-spheres and computes the
specification of the grid box for grid.
'''

def parse_line(line):
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

def min_and_max(tuple):
    '''
    Given a tuple (xcoord,ycoord,zcoord,radius) return another tuple
    with minimum and maximum coordinates (xmin,ymin,zmin,xmax,ymax,zmax).
    Xmin=xcoord-radius, xmax=xcoord+radius.
    '''
    (xcoord,ycoord,zcoord,radius)=tuple
    xmin=xcoord-radius
    xmax=xcoord+radius
    ymin=ycoord-radius
    ymax=ycoord+radius
    zmin=zcoord-radius
    zmax=zcoord+radius
    return (xmin,ymin,zmin,xmax,ymax,zmax)

def parse_file(lines):
    '''
    Go through all the lines in a file and return a list of 
    all 4-tuples.
    '''
    t4list=[]
    for line in lines:
        keyw=line[0:6]
        if keyw == b"ATOM  " or keyw == b"HETATM":
            tuple=parse_line(line)
            if tuple:
                t4list.append(tuple)
    return t4list

def gen_minmax(t4list):
    '''
    Go through the 4-tuple list and for each element 
    add an element to the 6-tuple list.
    '''
    t6list=[]
    for t4 in t4list:
        t6=min_and_max(t4)
        t6list.append(t6)
    return t6list

def find_minmax(t6list):
    '''
    Go through the list of 6-tuples and find the index for 
    the minimum and maximum coordinates. Return the result
    as a 12-tuple
    (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax).
    '''
    ixmin=0
    iymin=0
    izmin=0
    ixmax=0
    iymax=0
    izmax=0
    (xmin,ymin,zmin,xmax,ymax,zmax)=t6list[0]
    ielm=-1
    for t6 in t6list:
        ielm+=1
        (xmin1,ymin1,zmin1,xmax1,ymax1,zmax1)=t6
        if xmin1 < xmin:
            ixmin = ielm
            xmin  = xmin1
        if ymin1 < ymin:
            iymin = ielm
            ymin  = ymin1
        if zmin1 < zmin:
            izmin = ielm
            zmin  = zmin1
        if xmax1 > xmax:
            ixmax = ielm
            xmax  = xmax1
        if ymax1 > ymax:
            iymax = ielm
            ymax  = ymax1
        if zmax1 > zmax:
            izmax = ielm
            zmax  = zmax1
    return (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax)

def find_npts(t12,spacing,factor):
    '''
    Given the 12-tuple compute the number of points for each dimension.
    The number of points is computed from the length of a dimension
    and the grid spacing. The number of grid points must be a multiple of
    factor. Return the 3-tuple of numbers of points.
    '''
    (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax)=t12
    xlen=(xmax-xmin)/spacing
    ylen=(ymax-ymin)/spacing
    zlen=(zmax-zmin)/spacing
    xnpts=int(xlen/float(factor)+1.0)*int(factor)
    ynpts=int(ylen/float(factor)+1.0)*int(factor)
    znpts=int(zlen/float(factor)+1.0)*int(factor)
    return (xnpts,ynpts,znpts)
    
def find_center(t12):
    '''
    Given the 12-tuple compute the center of the receptor.
    Return the 3-tuple of center coordinates.
    '''
    (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax)=t12
    xcen=(xmax+xmin)/2.0
    ycen=(ymax+ymin)/2.0
    zcen=(zmax+zmin)/2.0
    return (xcen,ycen,zcen)
    
def find_lengths(t12):
    '''
    Given the 12-tuple compute the lengths of the receptor.
    Return the 3-tuple of lengths.
    '''
    (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax)=t12
    xlen=(xmax-xmin)
    ylen=(ymax-ymin)
    zlen=(zmax-zmin)
    return (xlen,ylen,zlen)

def parse_arguments():
    '''
    Parse the arguments.
    '''
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("filename",help="PQR file to generate grid-box from")
    parser.add_argument("--buffer",dest="buffer",help="Size of buffer region",default="0.0")
    args=parser.parse_args()
    return args

if __name__ == "__main__":
    args=parse_arguments()
    fobj=open(args.filename,"rb")
    buf=float(args.buffer)
    contents=fobj.readlines()
    fobj.close()
    t4list=parse_file(contents)
    t6list=gen_minmax(t4list)
    t12=find_minmax(t6list)
    npts=find_npts(t12,0.375,2)
    (xpts,ypts,zpts)=npts
    cent=find_center(t12)
    (xcen,ycen,zcen)=cent
    leng=find_lengths(t12)
    (xlen,ylen,zlen)=leng
    (ixmin,iymin,izmin,ixmax,iymax,izmax,xmin,ymin,zmin,xmax,ymax,zmax)=t12
    xlen+=2.0*buf
    ylen+=2.0*buf
    zlen+=2.0*buf
    xmin-=buf
    ymin-=buf
    zmin-=buf
    xmax+=buf
    ymax+=buf
    zmax+=buf
    print(f"HEADER    CORNERS OF BOX")
    print(f"REMARK    CENTER (X Y Z)   {xcen:.3f}  {ycen:.3f}  {zcen:.3f}")
    print(f"REMARK    DIMENSIONS (X Y Z)   {xlen:.3f}  {ylen:.3f}  {zlen:.3f}")
    print(f"ATOM      1  DUA BOX     1    {xmin:8.3f}{ymin:8.3f}{zmin:8.3f}")
    print(f"ATOM      2  DUB BOX     1    {xmax:8.3f}{ymin:8.3f}{zmin:8.3f}")
    print(f"ATOM      3  DUC BOX     1    {xmax:8.3f}{ymin:8.3f}{zmax:8.3f}")
    print(f"ATOM      4  DUD BOX     1    {xmin:8.3f}{ymin:8.3f}{zmax:8.3f}")
    print(f"ATOM      5  DUE BOX     1    {xmin:8.3f}{ymax:8.3f}{zmin:8.3f}")
    print(f"ATOM      6  DUF BOX     1    {xmax:8.3f}{ymax:8.3f}{zmin:8.3f}")
    print(f"ATOM      7  DUG BOX     1    {xmax:8.3f}{ymax:8.3f}{zmax:8.3f}")
    print(f"ATOM      8  DUH BOX     1    {xmin:8.3f}{ymax:8.3f}{zmax:8.3f}")
    print(f"CONECT    1    2    4    5")
    print(f"CONECT    2    1    3    6")
    print(f"CONECT    3    2    4    7")
    print(f"CONECT    4    1    3    8")
    print(f"CONECT    5    1    6    8")
    print(f"CONECT    6    2    5    7")
    print(f"CONECT    7    3    6    8")
    print(f"CONECT    8    4    5    7")
