# Dock-fix

The grid program in Dock fails to read OpenBabel's MOL2 file format
correctly. To fix this problem the protein MOL2 files were edited
by hand. This required 3 changes:

- Add one blank line after `GASTEIGER` and before `@<TRIPOS>ATOM`
- Increment the number of substructures from 0 to 1 on line 3
- Append a `@<TRIPOS>SUBSTRUCTURE` at the end of the file

Another issue is that Dock cannot handle substances that consist of
multiple disconnected fragments. This impacts drugs that are salts.
In most cases the code will seg-fault on such systems. In
`conf_gen_ag.cpp` a simple test was added to check that all
atoms were assigned to a segment. If not, then some atoms
could not be reached following the network of covalent bonds.
If the test fails then the code is now aborted. 

To install this fix copy `conf_gen_ag.cpp` to `dock6/src/dock`
and recompile.

The next issue is that Dock might try to align the ligand with 
the spheres that define the pocket. To do this the code tries to
allocate `squared((# spheres)x(# ligand atoms))`, but this dimension
is stored in an `int` (see line 738 in `orient.cpp`). This
can generate an overflow causing weird behavior. A test was added
to avoid the `int` overflow and abort the job with a sensible error
message.

To install this fix copy `orient.cpp` to `dock6/src/dock`
