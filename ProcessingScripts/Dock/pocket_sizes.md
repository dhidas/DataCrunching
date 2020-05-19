# Pocket sizes

This file contains the buffer dimensions needed around each pocket to generate
grids that are large enough for Dock6. The numbers in the table below are what is
needed in `example.sh` on the `prepare_target.sh` line.

Also note that `pqr2sph.py` now prunes the list of spheres to at most 100
spheres to avoid memory issues in Dock6.

Finally, in getting Dock6 to work I have also skipped the 100 largest molecules
from the list. These molecules were so large that Dock6 effectively cannot 
handle them. The Dock6 fix in Orient.cpp would probably help but only in as
far that it would generate a more informative error message.

| Pocket       | Buffer sizes |
| ------------ | ------------ |
| 3CLPro     1 | 14.0         |
| ADRP-ADPR  1 | 15.5         |
| ADRP-ADPR  5 | 18.0         |
| ADRP       1 | 15.0         |
| ADRP      12 | 16.0         |
| ADRP      13 | 16.0         |
| Nsp9       2 | 17.0         |
| Nsp9       7 | 17.0         |
| Nsp10      1 | 16.0         |
| Nsp10      3 | 14.0         |
| Nsp10     26 | 20.0         |
| Nsp15      1 | 17.0         |
| Nsp15      2 | 18.0         |
| PLPro      3 | 18.0         |
| PLPro      4 | 18.0         |
| PLPro      6 | 19.0         |
| PLPro     23 | 18.0         |
| PLPro     50 | 10.0         |
| CoV        1 | 18.0         |
| CoV        2 | 14.0         |
| Orf7a      2 | 18.0         |

