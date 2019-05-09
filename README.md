txt2mpig, genmpig, bfs-coco, ghs-coco
=====================================

`txt2mpig` converts the edge-centric graph model stored in a text file to a packed binary format.

`genmpig` generates an Erdos-Renyi graph and stores the edges in the same format as `txt2mpig`.

`bfs-coco` reads an input file generated by `txt2mpig` or `genmpig` and outputs the connected
components in the graph along with some basic component statistics. The connected components
are found by a breadth-first search (BFS) algorithm.

`ghs-coco` reads an input file generated by `txt2mpig` or `genmpig` and outputs the connected
components in the graph along with some basic component statistics. The connected components
are found by a modification of the Gallager-Humblet-Spira (GHS) algorithm for finding minimum
spanning trees.

txt2mpig
--------

### Usage

`mpirun -np` <node count> `txt2mpig -- -i` <input file> `-o` <output file> `-b` <edges per output block>

`-i` <input file>  : A textual edge-centric graph model representation with one edge per line.
                     Edge format: "[0-9]+ [0-9]+\n" (2 unsigned integers separated by whitespace)
`-o` <output file> : A binary edge-centric graph model representation with endpoint labels packed in tandem.
                     Edge format: <u|v> (4 bytes/label, 8 bytes/edge)
					 Default: `./graph.ecg`
`-b` <edges/block> : Output block buffering size in edges.
	                 Default: 1

### Notes
`txt2mpig` was designed to work with edge-centric text files available from the SNAP database
(https://snap.stanford.edu/data/). It reads the file in line by line, and discards lines that do not match
the expected format exactly (as understood by sscanf). The program does not deduplicate the input.

The edges are read at the root node until the read buffer is full or EOF is reached.
When this occurs, a portion of the read buffer is scattered to each participating machine which
then makes an independent write to a contiguous, non-overlapping area in the output file.

When the writes are complete, the absolute write offset is updated, and the read buffer offset is reset.
The transfer of data continues until the input file is exhausted. In the case where the read buffer
is not evenly divisible by the number of participating machines, the remainder of that division
is distributed to the lower ranked nodes. This remainder is accounted for in the per-machine
offset calculations so as to preserve the integrity of the output.

The default argument for the block buffering size is not a recommendation by any stretch.
Users should determine the best value for the particular use case. Best practice would seem to
be a multiple of the operating system block size divided by 8 bytes/edge.

genmpig
-------

### Usage

`mpirun -np` <node count> `genmpig -- -p` <edge probability> `-n` <population size> `-b` <edges/block>
                                      `-o` <output file> [`-m`]

`-o` <output file> : A binary edge-centric graph model representation with endpoint labels packed in tandem.
                     Edge format: <u|v> (4 bytes/label, 8 bytes/edge)
					 Default: `./graph.ecg`
`-b` <edges/block> : Output block buffering size in edges.
	                 Default: 1
`-p` <edge probability> : The independent probability of existence for any given edge in the graph.
                          The argument should be a double-precision value in the range [0,1].
                          Default: 0
`-n` <population size>  : The population size of the graph in terms of vertices.
	                      The argument should be a 32-bit unsigned integer.
                          Default: 1
`-m` : A switch to use Monte-Carlo edge selection instead of the naive Las Vegas edge selection.

### Notes
`genmpig` generates edge-centric graph model files in a packed binary format. A couple parameters are
the same or similar to those used by `txt2mpig`, but the write strategy is quite different since there
is an element of randomness introduced. In order to maximize the speed with which edges are generated,
each machine participates on behalf of a single node and its incident edges. To avoid duplicating effort,
each node under consideration only tests edges to nodes with higher labels since the edge to node
with lower label would have been generated (or not) during the consideration of the lower label node.

Additionally, there is no predictor that can calculate file offsets _a priori_ since each machine might
end up with an imbalanced load. Instead, a task-based parallel approach is employed and all machines
use a shared file pointer to accomplish writes via the MP/IO interface.

The `-p` and `-n` arguments are fairly self-explanatory in the default edge-selection context.
Under this scenario, each edge is tested by a Bernoulli sampling against a distribution with success
probability as given. This approach, while probabilistically sound, does not scale particularly well.
The combinatorial explosion of predicates becomes almost insurmountable.

To this end, the `-m` option is provided. Instead of doing individual Bernoulli tests, a random variable
with a binomial distribution is sampled to find the _number_ of incident edges against a given node.
The specific edges are taken from a random sampling of a uniform distribution of the graph's node labels.
Hence, whereas the naive selection scheme could be termed a Las Vegas randomized algorithm, this new scheme
falls under the Monte Carlo family of randomized algorithms.

The Monte Carlo edge selection scheme cannot guarantee that the random sampling of node labels will not
produce duplicate edges. For instance, if the edge probability is given as 1, the scheme may not generate
a complete graph. The naive scheme, on the other hand, will generate a complete graph.

Whether this massive increase in generation speed justifies the loss of probabilistic correctness remains
to be tested under experiment. It may be that with slight modifications, the probabilistic correctness
could be maintained while achieving much better efficiency than the naive approach.

Quick tests show that using the `-m` option usually generates a graph in subsecond timing or just
over a second, even up to `-n` = 2^20. The naive approach, conversely, can take up to almost 25 _minutes_
for a maximum tested `-n` = 2^17. For explicit details, see the genmpig*.log files in the project root.

NB: The program generates self-edges to ensure that singletons appear when loading the graph from the
file. Programs that read the packed binary format generated by `txt2mpig` and `genmpig` should handle
such pseudo-edges and deduplicate during the input phase.

bfs-coco
--------

### Usage

`mpirun -np` <node count> `bfs-coco -- -i`  <output file> [`-v`]

`-i` <input file> : A binary edge-centric graph model representation with endpoint labels packed in tandem.
                    Edge format: <u|v> (4 bytes/label, 8 bytes/edge)
                    Default: `./graph.ecg`
`-v` : Verbose mode. Dumps a list of graph nodes and their incident edges.

### Notes

`bfs-coco`, or *BFS*-*Co*nnected *Co*mponents, uses repeated BFS tree construction to search out connected
components in an arbitrary graph. The program loads the edge-centric binary input file and converts it in
memory into a vertex-centric format suitable for exploring via BFS.

Vertices are distributed among the participating machines according to a universal hash function described
in the COSC-6326 distributed algorithms class. The formula for the hashing is given as

	((_av_ + _b_) mod _p_) mod _n_

where _a_ is a random number in the range [1,_p_-1], _b_ is a random number in the range [0,_p_-1],
_n_ is the number of vertices in the graph, and _p_ is a prime number strictly greater than _n_.
For the sake of implementation, the Mersenne prime M61 (i.e., 2^61 - 1) is used in all cases for _p_.
The numbers _a_ and _b_ are randomly chosen at the start of each run of the program. The value for
_n_ is determined after the whole input file is read.

Since the input file size is known at runtime, and since the file format is embarrassingly parallel due
to the fixed size of each edge, all machines participate in reading a chunk of the file. The edges are
then distributed to the relevant machines according to the hash function described above. Once received
by the appropriate machine, the corresponding vertex is created at that machine. The convention in the
edge exchange messages is to place the target machine's node first. Again because of the fixed size,
there is no need to include headers in the messages since the size information is communicated via
`MPI_Gather` and the actual receipt of the edge information by `MPI_Gatherv`.

Once the vertex-centric model is constructed at each machine, a series of breadth-first searches
are executed until all vertices in the graph have been reached. The root of the current BFS tree, which
lends its name to the component, is chosen as the minimum of all ungrouped vertices across all machines.

The BFS root invites its neighbors to join the tree, who in turn invite their neighbors, and so on.
This is termed the broadcast phase.

Once a vertex who has joined a tree has heard a positive or negative response from all its neighbors,
it sends an upcast message itself to the parent vertex who first contacted it about joining a BFS tree.
The upcast message consists of the triple containing the parent ID, the child ID, and the number of
grandchildren (and further descendants) being contributed to the parent. In the case, where a node
receives an invitation to join after it had already accepted another invitation, it sends the same
upcast message with a count of zero (0). This is termed the upcast phase.

During this phase, nodes move from the input map to the output map. This is used as an efficient
termination condition since the algorithm terminates when the global input map is empty. Since
the program uses the C++ STL std::map under the hood, picking the minimum remaining vertex is
also efficient since the std::map iterator is ordered by key, which for implementation purposes
is the vertex ID.

The broadcast and upcast phases alternate synchronously across all machines until the BFS root
receives responses from all its children. At this point, the component labeled with the BFS root
is added to the component forest with the population accrued by the BFS root.

If the global input map is empty (i.e., all local input maps are empty), the algorithm terminates.
If not, the algorithm continues with the BFS root determination and another BFS search from that root.

The program went through several iterations to get to the current source. Frequent use of valgrind,
gdb, and callgrind were necessary to suss out the major inefficiencies, improper memory management,
and logic flaws that were encountered. In general, where memory could be sacrificed for the sake of
time, it was. The incoming edge map and the broadcast and upcast vertices vectors are prime examples.
These structures allowed for a limitation to be placed on node iteration.

Also, development went back and forth over destructive vs. non-destructive termination condition
checking. Ultimately, a non-destructive approach was employed for the sake of memory safety and for
providing a backbone which could be extended for other graph algorithms with the preservation of the
BFS tree in the output map. The preservation of the vertex neighbors set was probably the singlemost
helpful strategy in avoiding undercounting or overcounting the number of upcast responses received by
any given node. Inquisitive users may use `git log` or similar tools to observe the development history.

Graphs of the BFS experiments and other details not covered in this README are available in
the report `cosc6326-pa2-michael-yantosca.pdf`.

ghs-coco
--------

### Usage

`mpirun -np` <node count> `ghs-coco -- -i`  <output file> [`-v` `-vv` `-g`]

`-i` <input file> : A binary edge-centric graph model representation with endpoint labels packed in tandem.
                    Edge format: <u|v> (4 bytes/label, 8 bytes/edge)
                    Default: `./graph.ecg`
`-v`  : Verbose mode. Dumps a list of graph nodes and their incident edges.
`-vv` : Very verbose mode. Dumps a list of graph nodes and their incident edges plus messages sent and
        phase boundaries.
`-g`  : Debug mode. Triggers asserts if certain invariants are not honored. This was a holdover from
        development but was left in as an option for the interested user.

### Notes

`ghs-coco`, or *G*allager-*H*umblet-*S*pira-*Co*nnected *Co*mponents, uses the GHS algorithm for finding
minimum spanning trees (MSTs) to discover connected components in a component-breadth-first, vertex-depth-first
fashion. It has the benefit of working on multiple components simultaneously, but the current implementation
suffers from a unicast messaging model that would be better refactored into a broadcast model where possible.

The distribution of _n_ vertices among _k_ machines follows the implementation in `bfs-coco`. Once the model
is constructed, the program goes through a series of find and merge phases until there are no more joins to
be performed among the discovered components. After this series of phases is complete, the program performs
a component census with a gather centered at the machine with rank 0, which outputs the received results.

The post-hoc census is done to keep the code simple and maintainable. An online census would necessitate
many more messages because of the frequent reparenting of vertices during the merge phase. As implemented,
the census requires _O(n)_ additional messages for the downcast and convergecast in each component.

Some optimizations have been made to limit unnecessary communication. Since singleton components are
known by their lack of incident edges _ab initio_, they do not participate in the exchange phases but
are partitioned from the other components at the start. This addresses one of the primary failing of the
BFS algorithm implemented in `bfs-coco`. The same optimization has not been retrofitted to `bfs-coco`
due to time constraints.

To keep track of visited edges, the neighbor labels are separated into active and inactive sets.
Any minimum weight outgoing edge (MWOE) found during the find phase will necessarily
be moved from the active set to the inactive set. In the implementation used in the tests,
there is an unnecessary call to move edges involved in `new` messages into the inactive set,
but these edges will have already been moved on account of their use in `join` messages.
The excess lines were removed in commit `353437c3f6f919117a03ffc89a0df69ad09b721d`.
Quick smoke tests did not yield substantial differences in execution time, so the results plotted
in the report should be representative nonetheless.

As with `bfs-coco`, specifying the `-v` option on the command line will dump the actual whole graph
in human-readable format to `stderr`. By specifying the `-vv` option, the user will receive not only
the graph dump, but also all the messages exchanged by each machine during the course of the algorithm.
The user is recommended to specify the `mpirun` `--output-filename` option for best results in this case.

Another option which may be of interest to users is the `-g` or debug option. This option arms asserts
at critical points in the code to ensure that vertices are in expected states when sending and receiving
certain types of messages. There is a significant overhead incurred by using this feature, but it was
invaluable in tracing down various logic flaws during the course of development.

Graphs of the GHS experiments and other details not covered in this README are available in
the report `cosc6326-final-michael-yantosca.pdf`. Appendix A of the report contains graphs making
direct comparison of the empirical complexity between the BFS and GHS connected components algorithms.

Known Issues
------------

* `gen2mpig` `-m` option does not generate complete graphs when `-p` = 1.
* `gen2mpig` naive edge selection scheme has not been tested above `-n` = 2^17. Half an hour is long enough!
* `bfs-coco` does not perform well on disconnected graphs when scaling out to high values of _k_.
* `ghs-coco` does not perform well on denser graphs due to the strictly unicast messaging model.

Contact
-------

Please send comments, questions, and bug reports to Michael Yantosca via e-mail at
mike@archivarius.net.
