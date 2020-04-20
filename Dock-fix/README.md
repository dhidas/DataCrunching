# Dock-fix

The grid program in Dock fails to read OpenBabel's MOL2 file format
correctly. A fix for this problem has be included in the MOL2 reader
in `io_mol2.c`. 

Replace the file by the same name in `dock6/src/grid` with the
one provided here and rebuild the `grid` program.
