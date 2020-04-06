#!/bin/bash
#
# This script takes 1 argument:
#
# $1 - the target
# $2 - number of grid points
# $3 - grid center


[ -z "$1" ] && exit
[ -z "$2" ] && exit
[ -z "$3" ] && exit


#AREA=/sphenix/user/purschke/covid19
AREA=$PWD

module load autodock openbabel


LINE=$SLURM_ARRAY_TASK_ID
echo "$LINE"

declare -a fields

line=$(cat $AREA/ena+db.sorted | sed -n "${LINE}p")
echo "line = $line"
[ -z "$line" ] && exit

fields=($line)

smiles=${fields[0]}
id=${fields[1]}


TARGET=$1
NPTS=$2
GRIDCENTER=$3
NAME=$(basename $TARGET .pdb)
NAME="$NAME.pdbqt"

tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)-$id
mkdir -p $tmp_dir
mkdir  -p results
cp $NAME $tmp_dir
cd $tmp_dir


#pythonsh $AUTODOCKTOOLS_UTIL/prepare_receptor4.py -A checkhydrogens -r $TARGET -o $NAME

echo "Name = $NAME"

echo "$smiles" | obabel -h --gen3d -ismi -omol2 > $id.mol2
pythonsh $AUTODOCKTOOLS_UTIL/prepare_ligand4.py -l $id.mol2  -o $id.pdbqt
pythonsh $AUTODOCKTOOLS_UTIL/prepare_gpf4.py  -l $id.pdbqt -r $NAME -p npts=$NPTS -p gridcenter=$GRIDCENTER -o $id.gpf
pythonsh $AUTODOCKTOOLS_UTIL/prepare_dpf42.py -l $id.pdbqt -r $NAME -p ga_run=20 -o $id.dpf

autogrid4 -p $id.gpf -l $id.glg
autodock4 -p $id.dpf -l $id.dlg


mv $id.dlg $AREA/results/
cd ..
rm -rf $tmp_dir
