# Dock 6.9

Running a docking assay with Dock 6.9 proceeds in a number of steps. The
steps associated with setting up the receptor need to be done only once, 
whereas the docking steps need to be repeated for every ligand.

## Preparing the receptor

Preparing the receptor (protein-pocket combination) involves a few steps.

### 1. Converting the protein into the right format

For Dock the protein needs to be in the MOL2 format. We can use 
OpenBabel to convert the PDB file into a MOL2 file. 
```
obabel -h -ipdb $1.pdb -omol2 > $1.mol2
```
The grid box needs to be defined. In dock the grid box is specified in 
a PDB file that specifies the 8 corners of the box. Same as in Autodock
this information can be calculated from the grid box center and the side 
lengths.

The protein needs to be converted to the MOL2 format.

The ligand needs to be generated from a SMILES string.
The charges need to be calculated using Antechamber.
Save the molecule in MOL2 format.

The pockets are identified by spheres through a program called sphgen.
Produce a pdb or sph file with the selected spheres.

Making the docking grid.
Showbox generates the grid box.
Output ix rec_box.pdb listing the 8 corners of the grid box.
The generation of this information is essentially the same as in Autodock.
We can just take the center of grid box and the box lengths and construct
the rec_box.pdb file from that.
