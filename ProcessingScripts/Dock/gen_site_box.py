#!/usr/bin/env python3

def parse_arguments():
    '''
    Parse the arguments.
    '''
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("center",help="grid-box center \"xcenter,ycenter,zcenter\"")
    parser.add_argument("length",help="grid-box lengths \"xlength,ylength,zlength\"")
    args=parser.parse_args()
    return args

if __name__ == "__main__":
    args=parse_arguments()
    clist=args.center.split(",")
    llist=args.length.split(",")
    xcoord=float(clist[0])
    ycoord=float(clist[1])
    zcoord=float(clist[2])
    xlength=float(llist[0])
    ylength=float(llist[1])
    zlength=float(llist[2])
    xmax=xcoord+xlength/2.0
    xmin=xcoord-xlength/2.0
    ymax=ycoord+ylength/2.0
    ymin=ycoord-ylength/2.0
    zmax=zcoord+zlength/2.0
    zmin=zcoord-zlength/2.0
    print(f"HEADER    CORNERS OF BOX")
    print(f"REMARK    CENTER (X Y Z)   {xcoord:.3f}  {ycoord:.3f}  {zcoord:.3f}")
    print(f"REMARK    DIMENSIONS (X Y Z)   {xlength:.3f}  {ylength:.3f}  {zlength:.3f}")
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
