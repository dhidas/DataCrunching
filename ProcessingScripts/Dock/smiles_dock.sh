#!/bin/bash
#
# This script takes 4 argument:
#
# - $1.sph the spheres file (in sph format)
# - $2 the name of the file containing SMILES strings to process.
# ?- $3 the pocket center given as "xcoord,ycoord,zcoord"
# ?- $4 the pocket lengths given as "xlength,ylength,zlength"
#
#. /software/anaconda3/etc/profile.d/conda.sh
#conda activate py3
#
# Conversion of SMILES string to MOL2:
#
# - Done one at a time as OpenBabel might crash attempting this
#   and if that happens only 1 molecule is lost this way
#
if [ ${#0} -gt 14 ]
then
  # this command was called with an explicit path
  path=`readlink -f $0`
  path=`dirname $path`
else
  # this command is in the PATH
  path=`which smiles_dock.sh`
  path=`dirname $path`
fi
export MYCWD=`pwd`
declare -a fields
while IFS= read -r line
do
  fields=($line)
  smiles=${fields[0]}
  id=${fields[1]}
  if [ -f $id.dlg ]
  then
    continue
  fi
  mkdir $id
  cd $id
  $path/echo_smiles.py "$smiles" | obabel -h --gen3d --conformer --nconf 100 --score energy -ismi -omol2 > ${id}.mol2
  #$DOCK_HOME/bin/antechamber -i ${id}_obabel.mol2 -fi mol2 -o ${id}.mol2 -fo mol2 -c bcc -j 5 -s 2 
  #obabel -imol2 ${id}_ac.mol2 -omol2 > ${id}.mol2
  cat > anchor_and_grow.in <<EOF
conformer_search_type                                        flex
write_fragment_libraries                                     no
user_specified_anchor                                        no
limit_max_anchors                                            no
min_anchor_size                                              40
pruning_use_clustering                                       yes
pruning_max_orients                                          100
pruning_clustering_cutoff                                    100
pruning_conformer_score_cutoff                               25.0
pruning_conformer_score_scaling_factor                       1.0
use_clash_overlap                                            no
write_growth_tree                                            no
use_internal_energy                                          yes
internal_energy_rep_exp                                      12
internal_energy_cutoff                                       100.0
ligand_atom_file                                             $MYCWD/$id/$id.mol2
limit_max_ligands                                            no
skip_molecule                                                no
read_mol_solvation                                           no
calculate_rmsd                                               no
use_database_filter                                          no
orient_ligand                                                yes
automated_matching                                           yes
receptor_site_file                                           $MYCWD/$1.sph
max_orientations                                             500
critical_points                                              no
chemical_matching                                            no
use_ligand_spheres                                           no
bump_filter                                                  no
score_molecules                                              yes
contact_score_primary                                        no
contact_score_secondary                                      no
grid_score_primary                                           yes
grid_score_secondary                                         no
grid_score_rep_rad_scale                                     1
grid_score_vdw_scale                                         1
grid_score_es_scale                                          1
grid_score_grid_prefix                                       ../grid
multigrid_score_secondary                                    no
dock3.5_score_secondary                                      no
continuous_score_secondary                                   no
footprint_similarity_score_secondary                         no
pharmacophore_score_secondary                                no
descriptor_score_secondary                                   no
gbsa_zou_score_secondary                                     no
gbsa_hawkins_score_secondary                                 no
SASA_score_secondary                                         no
amber_score_secondary                                        no
minimize_ligand                                              yes
minimize_anchor                                              yes
minimize_flexible_growth                                     yes
use_advanced_simplex_parameters                              no
simplex_max_cycles                                           1
simplex_score_converge                                       0.1
simplex_cycle_converge                                       1.0
simplex_trans_step                                           1.0
simplex_rot_step                                             0.1
simplex_tors_step                                            10.0
simplex_anchor_max_iterations                                500
simplex_grow_max_iterations                                  500
simplex_grow_tors_premin_iterations                          0
simplex_random_seed                                          0
simplex_restraint_min                                        no
atom_model                                                   all
vdw_defn_file                                                $DOCK_HOME/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               $DOCK_HOME/parameters/flex.defn
flex_drive_file                                              $DOCK_HOME/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ../${id}_anchor_and_grow
write_orientations                                           no
num_scored_conformers                                        20
rank_ligands                                                 no
EOF
  dock6 -i anchor_and_grow.in -o ${id}_anchor_and_grow.out
  cd ..
  #rm -rf ${id}
done < $2
