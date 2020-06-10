#!/bin/bash
#
#  mlp 4/17 adapting to dock6
#  dah 4/20 adapt to NSLSII cluster
#
# This script takes 4 argument:
#
# $1 - the target mol2 file
# $2 - the PQR file
# $3 - the molecule number from ena-db, this is the line number
#      $2 is leading 0-padded, e.g. 000010 !  (formatted this way to ease the log file naming)
# $4 - nr of molecules (usually 10)
# $5 - the "stem" for filename generation e.g. PLPro_chainA
# $6 - the pocket number

# we use $5 and $6  to make a directory and file path, e.g.

# /sdcc/covid19/purschke/PLPro_chainA_Pocket1/PLPro_chainA_Pocket1.dlg
# (the name re-appears so it cannot be mistaken)
# also note that the directory needs to exist

echo "number of parameters : $#"

[ $# -lt 5 ] && exit
 
STARTTIME=$(date +%s)
echo  "starting at $STARTTIME"

echo $HOSTNAME

echo $*


MOL=$1
PQR=$2
#LINENR=$(echo $3 | sed 's/^0*//')
LINENR=$SLURM_ARRAY_TASK_ID
COUNT=$3
NAMESTEM=$4
POCKET=$5


AREA=/GPFS/APC/dhidas/DataCrunching/ProcessingScripts/BNL-CSI-NSLSII
OUTPUT=/GPFS/APC/dhidas/covid-19_data/dock6


#source /opt/sphenix/core/bin/sphenix_setup.sh -n
#source $AREA/setup.sh
module load gcc/7 mpich/3 dock python/3 openbabel
export PATH=/GPFS/APC/dhidas/DataCrunching/ProcessingScripts/Dock:$PATH

# let's make it so that we can test interactively
#[ -z "$_CONDOR_SCRATCH_DIR" ] && _CONDOR_SCRATCH_DIR=`pwd`/scratch
#[ -d "$_CONDOR_SCRATCH_DIR" ] || mkdir $_CONDOR_SCRATCH_DIR
#cd $_CONDOR_SCRATCH_DIR
_CONDOR_SCRATCH_DIR=$(mktemp -d -t ci-XXXXXXXXXX)
#_CONDOR_SCRATCH_DIR=/tmp/dhidas/a/
mkdir -p $_CONDOR_SCRATCH_DIR
cd $_CONDOR_SCRATCH_DIR

cp $MOL .
MOL=$(basename $MOL)

SPH=${PQR%.pqr}.sph

cp ${PQR%.pqr}_grid* .
GRID=${PQR%.pqr}_grid
GRID=$(basename $GRID)

cp $PQR .
PQR=$(basename $PQR)


MYCWD=$(pwd)

# we are generating (once) the target info 

pwd

############ done with general preparation


declare -a fields


TO=$(expr $LINENR + $COUNT - 1 ) 

#MNR = molecule number - we go (for now) through $COUNT
for MNR in $(seq $LINENR $TO) ; do 

    line=$(cat $AREA/ena+db.unique | sed -n "${MNR}p")
    echo "line nr $MNR - molecule = $line"  

    [ -z "$line" ] && exit

    fields=($line)

    smiles=${fields[0]}
    id=${fields[1]}

    echo "id = $id"
    # write this to sterr (separate file)  for easier failure analysis
    echo "line nr $MNR - molecule = $id"  1>&2

    mkdir  $_CONDOR_SCRATCH_DIR/$id
    cd  $_CONDOR_SCRATCH_DIR/$id

    cp $SPH .
    SPH=$(basename $SPH)
    echo "new SPH $SPH"
    pwd

    echo_smiles.py "$smiles" | obabel -h --gen3d --conformer --nconf 100 --score energy -ismi -omol2 > $id.mol2

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
receptor_site_file                                           $SPH
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
grid_score_grid_prefix                                       ../$GRID
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
ligand_outfile_prefix                                        ${id}_anchor_and_grow
write_orientations                                           no
num_scored_conformers                                        20
rank_ligands                                                 no
EOF

    dock6 -i anchor_and_grow.in -o ${id}_anchor_and_grow.out
    echo "dock done"

    obabel -l 1 -imol2 ${id}_anchor_and_grow_scored.mol2 -osdf | head --lines=-1 > $id.sdf
    d6_score=`grep "Grid_Score:" ${id}_anchor_and_grow.out | awk '{print $2}'`
    echo ">  <Dock>"  >> $id.sdf
    echo $d6_score    >> $id.sdf
    echo              >> $id.sdf
    echo ">  <TITLE>" >> $id.sdf
    echo $id          >> $id.sdf
    echo              >> $id.sdf
    echo "\$\$\$\$"   >> $id.sdf

    DIR="${OUTPUT}/${NAMESTEM}_Pocket${POCKET}"
    [ -d "DIR" ] || mkdir -p $DIR
    ls -l 

    cp ${id}.mol2 $DIR/
    cp $id.sdf $DIR/${id}.sdf
done

cd /tmp
/bin/rm -rf $_CONDOR_SCRATCH_DIR

ENDTIME=$(date +%s)
ET=$(date)
DURATION=$(expr $ENDTIME - $STARTTIME)
echo  "ending at $ENDTIME - $ET  Duration: $DURATION"

