# Dock 6.9

Running a docking assay with Dock 6.9 proceeds in a number of steps. The
steps associated with setting up the receptor need to be done only once, 
whereas the docking steps need to be repeated for every ligand.

Note that Dock 6.9 produces a MOL2 structure as an output. This file
contains the best docking pose and contains the corresponding docking
score in the comments.

## Preparing the receptor

Preparing the receptor (protein-pocket combination) involves a few steps
that can be performed with the `prepare_target.sh` script:
```
./prepare_target.sh <protein> <pocket>
```
Where `<protein>` is a PDB file (specified without the .pdb extension) 
containing the protein structure, and `<pocket>` is a PQR file
(specified without the .pqr extension) contain the list of sphere 
defining the pocket.

For example
```
./prepare_target.sh 3CLPro_protein pocket1
```

The steps involve the details given below.

### 1. Converting the protein into the right format

For Dock the protein needs to be in the MOL2 format. We can use 
OpenBabel to convert the PDB file into a MOL2 file. 
```
obabel -h -ipdb $1.pdb -omol2 > $1.mol2
```

### 2. Generate the pocket spheres file

The docking procedure itself Dock wants the spheres that locate the pocket.
Hence it makes sense to extract the piece of the PQR file corresponding to 
the pocket (manually) and then generate the sphere file that `dock` expects.
```
$path/pqr2sph.py "$2.pqr" > "$2.sph"
```

### 3. Generate the grid box file

The grid box needs to be defined. In dock the grid box is specified in 
a PDB file that specifies the 8 corners of the box. Same as in Autodock
this information can be calculated from the grid box center and the side 
lengths. However, for the docking procedure itself Dock also wants the
spheres that locate the pocket. Hence is makes sense to extract the piece
of the PQR file corresponding to the pocket and then calculate the grid 
box from there:
```
$path/gen_site_box.py "$2.pqr"  > site_box.pdb
```

### 4. Calculate the docking grids

The calculation of the grid is performed by the `grid` program of Dock.
Note that the grid calculation with Dock is a bit slower than the grid
calculation with Autodock.
```
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
box_file                       site_box.pdb
vdw_definition_file            $DOCK_HOME/parameters/vdw_AMBER_parm99.defn
score_grid_prefix              grid
EOF
$DOCK_HOME/bin/grid -i grid.in -o grid.out -v
```

## Dock a set of molecules

To dock a set of molecules two pieces of information are needed. These pieces
are the file with spheres as well as the file that contains the SMILES strings
of the molecules to be docked. The docking is done with the Dock version 
of the `smiles_dock.sh` script.
```
./smiles_dock.sh <pocket> <smiles>
```
where `<pocket>` is the file containing the spheres (specified without the .sph 
extension), the `<smiles>` file is the file containing lines with a SMILES strings
followed by an identifier. 

Example:
```
./smiles_dock.sh pocket1 ena+db-small.can
```

Executing this work involves a number of steps given below.

### 1. Convert SMILES string to MOL2

The SMILES string can be converted to a MOL2 file using OpenBabel
as was done with the Autodock approach.
```
echo "$smiles" | obabel -h --gen3d -ismi -omol2 > ${id}.mol2
```
The Dock tutorials recommend calculating the charges with the 
AM1-BCC approach rather then using the Gasteiger charges (that
OpenBabel would generate). However, Antechamber changes the atom
type names such that Dock won't recognize them while looking for
the Lennard-Jones potentials. Hence we stick with Openbabel.

### 2. Dock the molecule

Docking the molecule involves running Dock with a particular input file.
The `anchor_and_grow.in` input file contains many parameters. Below we
just give the ones with settings specific to our approach.
```
  cat > anchor_and_grow.in <<EOF
...
ligand_atom_file                                             $MYCWD/$id/$id.mol2
...
receptor_site_file                                           $MYCWD/$1.sph
...
grid_score_grid_prefix                                       ../grid
...
vdw_defn_file                                                $DOCK_HOME/parameters/vdw_AMBER_parm99.defn
flex_defn_file                                               $DOCK_HOME/parameters/flex.defn
flex_drive_file                                              $DOCK_HOME/parameters/flex_drive.tbl
ligand_outfile_prefix                                        ../${id}_anchor_and_grow
...
EOF
  dock6 -i anchor_and_grow.in -o ../${id}_anchor_and_grow.out
```

