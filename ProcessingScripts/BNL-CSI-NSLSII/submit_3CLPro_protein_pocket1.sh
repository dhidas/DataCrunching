#!/bin/bash

TARGET=/GPFS/APC/dhidas/DataCrunching/Data/TARGETS/3CLPro_pockets/3CLPro_protein.mol2
PQR=/GPFS/APC/dhidas/DataCrunching/Data/TARGETS/3CLPro_pockets/pocket1_vert.pqr
NPTS=1
NAME=3CLPro
POCKET=1

LDIR=/GPFS/APC/dhidas/dock6_logs/${NAME}_Pocket${POCKET}
mkdir -pv ${LDIR}

#1-310693
sbatch -p covid -J ${NAME}${POCKET} -a 1-310693 -o ${LDIR}/${NAME}_Pocket${POCKET}_%a.out -e ${LDIR}/${NAME}_Pocket${POCKET}_%a.err --wrap="./run_dock6.sh ${TARGET} ${PQR} ${NPTS} ${NAME} ${POCKET}"
