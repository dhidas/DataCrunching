#!/bin/bash
#
# Set the target up for Dock 6.9
#
# This script takes 3 argument:
#
# - $1.mol2 is the file containing the structure of the protein
# - $2.pqr  is the file containing the fpocket alpha-spheres
# - $3      is the size of the buffer around the pocket for
#           defining the grids
#
if [ ${#0} -gt 17 ]
then
  # this command was called with an explicit path
  path=`dirname $0`
else
  # this command is in the PATH
  path=`which prepare_target.sh`
  path=`dirname $path`
fi
buf=0.0
if [ $# -eq 3 ]
then
  buf=$3
fi
export MYCWD=`pwd`
#
# OpenBabel generates MOL2 files the Dock6's grid program cannot read.
#
#obabel -h -ipdb "$1.pdb" -omol2 | sed 's/GASTEIGER/GASTEIGER\n/' > "$1.mol2"
#cat <<EOF >> "$1.mol2"
#@<TRIPOS>SUBSTRUCTURE
#1 molecule PERM 0 **** **** 0 ROOT
#EOF
$path/pqr2sph.py "$2.pqr" > "$2.sph"
$path/gen_site_box.py "$2.pqr" --buffer "$buf" > site_box.pdb
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
receptor_file                  $1.mol2
box_file                       $MYCWD/site_box.pdb
vdw_definition_file            $DOCK_HOME/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              grid
EOF
$DOCK_HOME/bin/grid -i grid.in -o grid.out -v
