#!/usr/bin/env python3

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

def process_lines(list_in):
    '''
    Go through all lines and convert each one to the spheres format, add it to the result,
    and return the list of output lines.
    '''
    list_out=[]
    length=len(list_in)
    list_out.append(f"DOCK spheres within 10.0 ang of ligands")
    list_out.append(f"cluster     1   number of spheres in cluster  {length:4d}")
    for line_in in list_in:
        keyw=line_in[0:6]
        if keyw == b"ATOM  " or keyw == b"HETATM":
            data=parse_pqr_line(line_in)
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
