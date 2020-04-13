#!/bin/bash
#
# Set the target up for Dock 6.9
#
# This script takes 3 argument:
#
# - $1.pdb is the file containing the structure of the protein
# - $2     is the grid box center  "xcoord,ycoord,zcoord"
# - $3     is the grid box lengths "xlength,ylength,zlength"
#
export MYCWD=`pwd`
obabel -h -ipdb $1.pdb -omol2 > $1.mol2
gen_site_box.py "$2" "$3" > site_box.pdb
cat > grid.in <<EOF
compute_grids                  yes
grid_spacing                   0.3
output_molecule                no
contact_score                  no
energy_score                   yes
energy_cutoff_distance         9999
atom_model                     a
attractive_exponent            6
repulsive_exponent             12
distance_dielectric            yes
dielectric_factor              4
bump_filter                    yes
bump_overlap                   0.75
receptor_file                  $MYCWD/$1.mol2
box_file                       $MYCWD/site_box.pdb
vdw_definition_file            $DOCK_PREFIX/dock6/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              grid
EOF
$DOCK_PREFIX/dock6/bin/grid -i grid.in -v
