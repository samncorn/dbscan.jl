# dbscan.jl

[![Build Status](https://github.com/samncorn/dbscan.jl/workflows/CI/badge.svg)](https://github.com/samncorn/dbscan.jl/actions)

Implementation of the algorithm from M. M. A. Patwary, D. Palsetia, A. Agrawal, W. -k. Liao, F. Manne and A. Choudhary, 
"A new scalable parallel DBSCAN algorithm using the disjoint-set data structure," SC '12: Proceedings of the 
International Conference on High Performance Computing, Networking, Storage and Analysis, Salt Lake City, UT, USA, 2012, 
pp. 1-11, doi: 10.1109/SC.2012.9.

Uses a KDTree for neighborhood queries.

Multithreading is planned but not yet fully implemented.

See /demos/ for an example of use.