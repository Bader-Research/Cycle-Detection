# Cycle-Detection: A Parallel Algorithm for Cycle Detection in Planar Partitioned Digraphs

This package provides an MPI implementation of a new, parallel
algorithm for detecting cycles in partitioned, planar directed graphs
that is both scalable in the graph and machine size, and performs well
in practice. As an example, on a p = 64 processor cluster, we have
solved an extremely large and difficult input problem with n = 228
vertices in less than five minutes.

References:

D. A. Bader, "A Practical Parallel Algorithm for Cycle Detection in
Partitioned Digraphs," September 1999. (Technical Report
UNM/AHPCC99-013).
