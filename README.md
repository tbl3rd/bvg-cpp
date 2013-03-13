bvg-cpp
=======

Read a mutation probability and gene data source file name from the
command line.  Initialize a bitvector population from the data
stream, and calculate the relation distance between all pairs of
vectors in the population.

Find an undirected graph that spans all the bitvectors and
minimizes the bit difference between adjacent bitvectors normalized
to the number of expected mutations.

Orient the graph into a rooted tree by first transforming the
spanning graph's vector and edge representation into a neighborhood
representation.  The leaves of the tree are the bitvectors with
only a single neighbor.  Find the leaves, note their parents, then
trim them from the graph, and repeat until there is a single
bitvector left, the "progenitor".

Run 'make test' to build the executables, run the tests, and check
the results.  This code was developed on MacOSX, but should run on
any standardly-endowed Unix system with a Makefile tweak or two.

Or just run 'make bitvectors-parents.data' for the "solution".
