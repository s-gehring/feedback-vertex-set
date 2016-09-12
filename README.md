# Compiling

Use the makefile. 

    make          


# Program usage

    ./prog < [path_to_graph_file]

With a graph file being a list of string tuples corresponding to edges, each tuple terminated by a line break and each element of said tuples separated by a white-space.
Comment lines begin with a '#' literal.


# Documentation 

1. General Case

First, we apply basic preprocessing procedures as deleting all vertices with self-loops and
deleting all vertices with degree strictly smaller than two. The next step is to iteratively
remove semidisjoint cycles, which are cycles from which at most one vertex has degree
greater two. There exists always an optimal solution containing this higher degree vertex.
Afterwards the graph is cut along all its containing bridges, potentially creating new
connected components. Notice that a bridge is never part of a cycle and thus is safe to
remove. We solve the feedback-vertex-set for each connected component separately and
union all local solutions. After this preprocessing we check whether the degree of each
vertex of the graph is three or smaller since these instances can be solved in polynomial
time by applying the procedure explained below. If the resulting graph contains vertices
with a higher degree, we furthermore try to minimize it via edge contractions. For this
step we take a vertex with degree two and replace it by an edge connecting its neighbors,
possibly creating multiedges. Notice that for every such multiedge, we need at least one
of its vertices. Since some multiedges may be overlapping, the last preprocessing step is
to create a graph consisting only of the multiedges and detect all minimal vertex covers
by a recursive procedure. These minimal vertex covers correspond to all the dierent
multiedge branchings which need to be considered and all of them are tested for nding
the best solution.
Now, the main part of the algorithm is applied: the iterative compression algorithm
which is presented in [3]. Within this algorithm, the approximative solution is computed
with the algorithm given in [1]. There is only one minor change: Instead of taking
the whole semidisjoint cycle into the unreduced approximate solution in every step, we
only take the vertex which has degree greater than two instead of taking the whole
semidisjoint cycle. It should be obvious that this is still optimal. Whenever we make
use of algorithm 1 given in [3], we first solve all the vertices with highest degree until
there are only vertices with degree smaller or equal three left. The reason for this is that
we want to apply the polynomial time degree three case algorithm as often as possible.

2. Degree Three Case

We detect and solve the degree three special case in polynomial time as described in [2].
We call an instance a degree three case, if in the problem DISJOINT-􀀀FVS(G; V1; V2; k)
all vertices of V1 have degree three or smaller, where V2 is a FVS of G of size k + 1 and
we want to compute a FVS completely in V1 of size k. The DISJOINT-FVS problem
appears in the iterative compression.
The key idea is the following: If we consider a G[V2]-spanning tree T, every edge in
G-E[T] creates a cycle. A V1 􀀀 adjacency 􀀀 matching is a partition of those edge in
G􀀀-E[T] in groups of one or two edges and every group of two edges share an endpoint
in V1. We are searching for a G[V2]-spanning tree that maximizes the number of two-
groups. Then we put for every two-group the shared vertex and for every one-group
an arbitrary vertex into the FVS. This is showed to be an optimal FVS in [2]. To get
this best G[V2]-spanning tree T we reduce this problem to a cographic matroid parity
problem of an auxiliary graph described in Section 3 of the paper. The cographic version
is itself solved by a linear matroid parity algorithm as described in [4]. To get a huge
speedup in the algorithm, we use the small rank update formula by Sherman, Morrison
and Woodbury to check if a matrix is non-singular. In every call of the degree three
routine we plug in random numbers of a finite field with size 2^64 for the indeterminates
and calculate the maximum full rank submatrix. If the full rank submatrix is not the
entire matrix we generate new random numbers and check if we get the same result. If
not, we repeat the whole routine. The size of the matrix is the difference of edges to the
nodes, this is smaller than the number of edges. We know that the instances have at
most 90000 edges, thus the size of the matrix is also bounded by this number. So by the
Schwarz-Zippel Lemma the failure probability is bounded by 90000/(2^64). Since we generate
random numbers twice in every call, we can bound it to (90000/(2^64))^2. The clock signal of
the processor is 3.6 GHz, therefore it can perform 3.6*10^9*1800<10^15 operations in 30
minutes, hence our routine is executed less then 1015 times in 30 minutes. The overall
failure probability is bounded by

    1-(1-(90000/2^64)^2)^10^15<10^-12
    
3. References

[1] V. Bafna, P. Berman, and T. Fujito. A 2-approximation algorithm for the undirected
feedback vertex set problem. SIAM Journal on Discrete Mathematics, 12(3):289-297,
1999.
[2] Y. Cao, J. Chen, and Y. Liu. On feedback vertex set, new measure and new struc-
tures. CoRR, abs/1004.1672, 2010.
[3] J. Chen, F. V. Fomin, Y. Liu, S. Lu, and Y. Villanger. Improved algorithms for
feedback vertex set problems. Journal of Computer and System Sciences, 74(7):1188-
1198, 2008.
[4] H. Y. Cheung, L. C. Lau, and K. M. Leung. Algebraic algorithms for linear matroid
parity problems. ACM Trans. Algorithms, 10(3):10:1-10:26, May 2014.