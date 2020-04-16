#!/bin/bash

./prepare_target.sh 3CLPro_protein pocket1
./smiles_dock.sh pocket1 ena+db-small.can
./summarize.sh | tee summary.txt
