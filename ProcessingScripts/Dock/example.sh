#!/bin/bash

rm -rf DB-2133 DB-2858 DB-3969 DB-7596 DB-9218
rm -f  grid.nrg grid.bmp grid.in grid.out
rm -f  DB-2133_anchor_and_grow_conformers.mol2 DB-2133_anchor_and_grow_scored.mol2
rm -f  DB-2858_anchor_and_grow_conformers.mol2 DB-2858_anchor_and_grow_scored.mol2
rm -f  DB-3969_anchor_and_grow_conformers.mol2 DB-3969_anchor_and_grow_scored.mol2
rm -f  DB-7596_anchor_and_grow_conformers.mol2 DB-7596_anchor_and_grow_scored.mol2
rm -f  DB-9218_anchor_and_grow_conformers.mol2 DB-9218_anchor_and_grow_scored.mol2
rm -f  summary.txt

./prepare_target.sh 3CLPro_protein pocket1
./smiles_dock.sh pocket1 ena+db-small.can
./summarize.sh | tee summary.txt
