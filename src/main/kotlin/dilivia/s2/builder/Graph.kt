package dilivia.s2.builder

import dilivia.s2.Assertions
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.Assertions.assertLT
import dilivia.s2.Assertions.assertTrue
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.assign
import dilivia.s2.builder.S2Builder.IsFullPolygonPredicate
import dilivia.s2.remove

// Identifies a vertex in the graph.  Vertices are numbered sequentially
// starting from zero.
typealias VertexId = Int

// A loop consisting of a sequence of edges.
typealias EdgeLoop = List<EdgeId>

typealias DirectedComponent = List<EdgeLoop>
typealias UndirectedComponent = Pair<List<EdgeLoop>, List<EdgeLoop>>
typealias EdgePolyline = List<EdgeId>

// A struct for sorting the incoming and outgoing edges around a vertex "v0".
data class VertexEdge(
    val incoming: Boolean,       // Is this an incoming edge to "v0"?
    val index: EdgeId,           // Index of this edge in "edges_" or "in_edge_ids"
    val endpoint: VertexId,      // The other (not "v0") endpoint of this edge
    val rank: Int                // Secondary key for edges with the same endpoint
)

// An S2Builder::Graph represents a collection of snapped edges that is passed
// to a Layer for assembly.  (Example layers include polygons, polylines, and
// polygon meshes.)  The Graph object does not own any of its underlying data;
// it is simply a view of data that is stored elsewhere.  You will only
// need this interface if you want to implement a new Layer subtype.
//
// The graph consists of vertices and directed edges.  Vertices are numbered
// sequentially starting from zero.  An edge is represented as a pair of
// vertex ids.  The edges are sorted in lexicographic order, therefore all of
// the outgoing edges from a particular vertex form a contiguous range.
//
// S2Builder::Graph is movable and copyable.  Note that although this class
// does not own the underlying vertex and edge data, S2Builder guarantees that
// all Graph objects passed to S2Builder::Layer::Build() methods will remain
// valid until all layers have been built.
//
// TODO(ericv): Consider pulling out the methods that are helper functions for
// Layer implementations (such as GetDirectedLoops) into s2builderutil_graph.h.
// Note that most of the parameters are passed by const reference and must
// exist for the duration of the Graph object.  Notes on parameters:


// "label_set_lexicon":
//   - a class that maps a LabelSetId to a set of S2Builder::Labels.
// "is_full_polygon_predicate":
//   - a predicate called to determine whether a graph consisting only of
//     polygon degeneracies represents the empty polygon or the full polygon
//     (see s2builder.h for details).
class Graph(
        // "options":
        //    - the GraphOptions used to build the Graph.  In some cases these
        //      can be different than the options provided by the Layer.
        val options: GraphOptions = GraphOptions(),

        // "vertices":
        //   - a vector of S2Points indexed by VertexId.
        val vertices: List<S2Point> = emptyList(),

        // "edges":
        //   - a vector of VertexId pairs (sorted in lexicographic order)
        //     indexed by EdgeId.
        val edges: List<Edge> = emptyList(),

        // "input_edge_id_set_ids":
        //   - a vector indexed by EdgeId that allows access to the set of
        //     InputEdgeIds that were mapped to the given edge, by looking up the
        //     returned value (an InputEdgeIdSetId) in "input_edge_id_set_lexicon".
        val inputEdgeIdSetIds: List<InputEdgeIdSetId> = emptyList(),

        // "input_edge_id_set_lexicon":
        //   - a class that maps an InputEdgeIdSetId to a set of InputEdgeIds.
        val inputEdgeIdSetLexicon: IdSetLexicon,

        // "label_set_ids":
        //   - a vector indexed by InputEdgeId that allows access to the set of
        //     labels that were attached to the given input edge, by looking up the
        //     returned value (a LabelSetId) in the "label_set_lexicon".
        val labelSetIds: List<LabelSetId> = emptyList(),

        val labelSetLexicon: IdSetLexicon,

        val is_full_polygon_predicate: IsFullPolygonPredicate

) {

    private var numVertices: Int = -1;  // Cached to avoid division by 24.

    // Returns the number of vertices in the graph.
    fun numVertices(): VertexId = numVertices

    // Returns the vertex at the given index.
    fun vertex(v: VertexId): S2Point = vertices[v]

    // Returns the total number of edges in the graph.
    fun numEdges(): EdgeId = edges.size

    // Returns the endpoints of the given edge (as vertex indices).
    fun edge(e: EdgeId): Edge = edges[e]

    // Returns a vector of edge ids sorted in lexicographic order by
    // (destination, origin).  All of the incoming edges to each vertex form a
    // contiguous subrange of this ordering.
    fun getInEdgeIds(): MutableList<EdgeId> {
        val inEdgeIds = mutableListOf<EdgeId>()
        repeat(numEdges()) { i -> inEdgeIds.add(i) }
        inEdgeIds.sortWith { ai, bi -> if (stableLessThan(reverse(edge(ai)), reverse(edge(bi)), ai, bi)) -1 else 1 }
        return inEdgeIds
    }

    // Given a graph such that every directed edge has a sibling, returns a map
    // from EdgeId to the sibling EdgeId.  This method is identical to
    // GetInEdgeIds() except that (1) it requires edges to have siblings, and
    // (2) undirected degenerate edges are grouped together in pairs such that
    // one edge is the sibling of the other.  Handles duplicate edges correctly
    // and is also consistent with GetLeftTurnMap().
    //
    // REQUIRES: An option is chosen that guarantees sibling pairs:
    //     (options.sibling_pairs() == { REQUIRE, CREATE } ||
    //      options.edge_type() == UNDIRECTED)
    fun getSiblingMap(): List<EdgeId> {
        val inEdgeIds = getInEdgeIds()
        makeSiblingMap(inEdgeIds)
        return inEdgeIds
    }

    // Like GetSiblingMap(), but constructs the map starting from the vector of
    // incoming edge ids returned by GetInEdgeIds().  (This operation is a no-op
    // except unless undirected degenerate edges are present, in which case such
    // edges are grouped together in pairs to satisfy the requirement that every
    // edge must have a sibling edge.)
    fun makeSiblingMap(in_edge_ids: MutableList<EdgeId>): Unit {
        Assertions.assert {
            (options.sibling_pairs == SiblingPairs.REQUIRE ||
                    options.sibling_pairs == SiblingPairs.CREATE ||
                    options.edge_type == S2Builder.EdgeType.UNDIRECTED)
        }
        repeat(numEdges()) { e -> Assertions.assert { (edge(e) == reverse(edge(in_edge_ids[e]))) } }
        if (options.edge_type == S2Builder.EdgeType.DIRECTED) return
        if (options.degenerate_edges == DegenerateEdges.DISCARD) return

        var e = 0
        while (e < numEdges()) {
            val v = edge(e).first
            if (edge(e).second == v) {
                assertLT(e + 1, numEdges())
                assertEQ(edge(e + 1).first, v)
                assertEQ(edge(e + 1).second, v)
                assertEQ(in_edge_ids[e], e)
                assertEQ(in_edge_ids[e + 1], e + 1)
                in_edge_ids[e] = e + 1
                in_edge_ids[e + 1] = e
                ++e
            }
            ++e
        }
    }

    // Returns the set of input edge ids that were snapped to the given
    // edge.  ("Input edge ids" are assigned to input edges sequentially in
    // the order they are added to the builder.)  For example, if input
    // edges 2 and 17 were snapped to edge 12, then input_edge_ids(12)
    // returns a set containing the numbers 2 and 17.  Example usage:
    //
    //   for (InputEdgeId input_edge_id : g.input_edge_ids(e)) { ... }
    //
    // Please note the following:
    //
    //  - When edge chains are simplified, the simplified edge is assigned all
    //    the input edge ids associated with edges of the chain.
    //
    //  - Edges can also have multiple input edge ids due to edge merging
    //    (if DuplicateEdges::MERGE is specified).
    //
    //  - Siblings edges automatically created by EdgeType::UNDIRECTED or
    //    SiblingPairs::CREATE have an empty set of input edge ids.  (However
    //    you can use a LabelFetcher to retrieve the set of labels associated
    //    with both edges of a given sibling pair.)
    fun inputEdgeIds(e: EdgeId): IdSetLexicon.IdSet = inputEdgeIdSetLexicon.idSet(inputEdgeIdSetIds[e])

    // Low-level method that returns an integer representing the entire set of
    // input edge ids that were snapped to the given edge.  The elements of the
    // IdSet can be accessed using input_edge_id_set_lexicon().
    fun inputEdgeIdSetId(e: EdgeId): InputEdgeIdSetId = inputEdgeIdSetIds[e]

    // Returns the minimum input edge id that was snapped to this edge, or -1 if
    // no input edges were snapped (see SiblingPairs::CREATE).  This is
    // useful for layers that wish to preserve the input edge ordering as much
    // as possible (e.g., to ensure idempotency).
    fun minInputEdgeId(e: EdgeId): InputEdgeId {
        val idSet = inputEdgeIds(e)
        return if (idSet.size() == 0) kNoInputEdgeId else idSet.values.first()
    }

    // Returns a vector containing the minimum input edge id for every edge.
    // If an edge has no input ids, kNoInputEdgeId is used.
    fun getMinInputEdgeIds(): List<InputEdgeId> {
        val minInputIds = mutableListOf<InputEdgeId>()
        repeat(numEdges()) { e -> minInputIds[e] = minInputEdgeId(e) }
        return minInputIds
    }

    // Returns a vector of EdgeIds sorted by minimum input edge id.  This is an
    // approximation of the input edge ordering.x
    fun getInputEdgeOrder(min_input_edge_ids: List<InputEdgeId>): List<EdgeId> {
        val order = mutableListOf<EdgeId>()
        min_input_edge_ids.indices.forEach { i -> order.add(i) }
        order.sortWith { a: EdgeId, b: EdgeId ->
            if (min_input_edge_ids[a] == min_input_edge_ids[b]) a.compareTo(b)
            else min_input_edge_ids[a].compareTo(min_input_edge_ids[b])
        }
        return order;
    }

    // Returns the set of labels associated with a given input edge.  Example:
    //   for (Label label : g.labels(input_edge_id)) { ... }
    fun labels(e: InputEdgeId): IdSetLexicon.IdSet = labelSetLexicon.idSet(labelSetIds[e])

    // Low-level method that returns an integer representing the set of
    // labels associated with a given input edge.  The elements of
    // the IdSet can be accessed using label_set_lexicon().
    fun labelSetId(e: InputEdgeId): LabelSetId = labelSetIds[e]

    // Convenience method that calls is_full_polygon_predicate() to determine
    // whether a graph that consists only of polygon degeneracies represents the
    // empty polygon or the full polygon (see s2builder.h for details).
    fun isFullPolygon(error: S2Error): Boolean = is_full_polygon_predicate.test(this, error)

    // Builds loops from a set of directed edges, turning left at each vertex
    // until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
    // (for LoopType::CIRCUIT) is found.  (Use LoopType::SIMPLE if you intend to
    // construct an S2Loop.)
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // If any degenerate edges are present, then each such edge is treated as a
    // separate loop.  This is mainly useful in conjunction with
    // options.degenerate_edges() == DISCARD_EXCESS, in order to build polygons
    // that preserve degenerate geometry.
    //
    // REQUIRES: options.degenerate_edges() == {DISCARD, DISCARD_EXCESS}
    // REQUIRES: options.edge_type() == DIRECTED
    fun getDirectedLoops(loop_type: LoopType, loops: MutableList<EdgeLoop>, error: S2Error): Boolean {
        assertTrue(options.degenerate_edges == DegenerateEdges.DISCARD ||
                options.degenerate_edges == DegenerateEdges.DISCARD_EXCESS)
        assertTrue(options.edge_type == S2Builder.EdgeType.DIRECTED)

        val leftTurnMap = mutableListOf<EdgeId>()
        if (!getLeftTurnMap(getInEdgeIds(), leftTurnMap, error)) return false
        val minInputIds = getMinInputEdgeIds()

        // If we are breaking loops at repeated vertices, we maintain a map from
        // VertexId to its position in "path".
        val path_index = mutableListOf<Int>()
        if (loop_type == LoopType.SIMPLE) repeat(numVertices) { path_index.add(-1) }

        // Visit edges in arbitrary order, and try to build a loop from each edge.
        val path = mutableListOf<EdgeId>()
        for (start in 0 until numEdges()) {
            if (leftTurnMap[start] < 0) continue

            // Build a loop by making left turns at each vertex until we return to
            // "start".  We use "left_turn_map" to keep track of which edges have
            // already been visited by setting its entries to -1 as we go along.  If
            // we are building vertex cycles, then whenever we encounter a vertex that
            // is already part of the path, we "peel off" a loop by removing those
            // edges from the path so far.
            var e = start
            var next: EdgeId
            while (leftTurnMap[e] >= 0) {
                path.add(e)
                next = leftTurnMap[e];
                leftTurnMap[e] = -1;
                if (loop_type == LoopType.SIMPLE) {
                    path_index[edge(e).first] = path.size - 1
                    val loop_start = path_index[edge(e).second]
                    if (loop_start < 0) continue
                    // Peel off a loop from the path.
                    val loop = path.subList(loop_start, path.size)
                    while (path.size > loop_start) path.removeLast()
                    path.remove(loop_start, path.size)
                    for (e2 in loop) path_index[edge(e2).first] = -1
                    canonicalizeLoopOrder(minInputIds, loop)
                    loops.add(loop)
                }
                e = next
            }
            if (loop_type == LoopType.SIMPLE) {
                assertTrue(path.isEmpty())  // Invariant.
            } else {
                canonicalizeLoopOrder(minInputIds, path)
                loops.add(path)
                path.clear()
            }
        }
        canonicalizeVectorOrder(minInputIds, loops)
        return true
    }

    // Returns a map "m" that maps each edge e=(v0,v1) to the following outgoing
    // edge around "v1" in clockwise order.  (This corresponds to making a "left
    // turn" at the vertex.)  By starting at a given edge and making only left
    // turns, you can construct a loop whose interior does not contain any edges
    // in the same connected component.
    //
    // If the incoming and outgoing edges around a vertex do not alternate
    // perfectly (e.g., there are two incoming edges in a row), then adjacent
    // (incoming, outgoing) pairs are repeatedly matched and removed.  This is
    // similar to finding matching parentheses in a string such as "(()())()".
    //
    // For sibling edge pairs, the incoming edge is assumed to immediately
    // follow the outgoing edge in clockwise order.  Thus a left turn is made
    // from an edge to its sibling only if there are no other outgoing edges.
    // With respect to the parentheses analogy, a sibling pair is ")(".
    // Similarly, if there are multiple copies of a sibling edge pair then the
    // duplicate incoming and outgoing edges are sorted in alternating order
    // (e.g., ")()(").
    //
    // Degenerate edges (edges from a vertex to itself) are treated as loops
    // consisting of a single edge.  This avoids the problem of deciding the
    // connectivity and ordering of such edges when they share a vertex with
    // other edges (possibly including other degenerate edges).
    //
    // If it is not possible to make a left turn from every input edge, this
    // method returns false and sets "error" appropriately.  In this situation
    // the left turn map is still valid except that any incoming edge where it
    // is not possible to make a left turn will have its entry set to -1.
    //
    // "in_edge_ids" should be equal to GetInEdgeIds() or GetSiblingMap().
    fun getLeftTurnMap(in_edge_ids: List<EdgeId>, left_turn_map: MutableList<EdgeId>, error: S2Error): Boolean {
        left_turn_map.assign(numEdges(), -1)
        if (numEdges() == 0) return true

        // Declare vectors outside the loop to avoid reallocating them each time.
        val v0_edges = mutableListOf<VertexEdge>()
        val e_in = mutableListOf<EdgeId>()
        val e_out = mutableListOf<EdgeId>()

        // Walk through the two sorted arrays of edges (outgoing and incoming) and
        // gather all the edges incident to each vertex.  Then we sort those edges
        // and add an entry to the left turn map from each incoming edge to the
        // immediately following outgoing edge in clockwise order.
        var out = 0
        var input = 0
        var out_edge = edge(out)
        var in_edge = edge(in_edge_ids[input])
        val sentinel = Edge(numVertices(), numVertices())
        var min_edge = min(out_edge, reverse(in_edge))
        while (min_edge != sentinel) {
            // Gather all incoming and outgoing edges around vertex "v0".
            val v0: VertexId = min_edge.first
            while (min_edge.first == v0) {
                val v1: VertexId = min_edge.second
                // Count the number of copies of "min_edge" in each direction.
                val out_begin = out
                var in_begin = input
                while (out_edge == min_edge) {
                out_edge = if(++out == numEdges()) sentinel else edge(out)
            }
                while (reverse(in_edge) == min_edge) {
                    in_edge = if(++input == numEdges()) sentinel else edge(in_edge_ids[input])
                }
                if (v0 != v1) {
                    addVertexEdges(out_begin, out, in_begin, input, v1, v0_edges)
                } else {
                    // Each degenerate edge becomes its own loop.
                    while (in_begin < input) {
                        left_turn_map[in_begin] = in_begin
                        ++in_begin
                    }
                }

                min_edge = min(out_edge, reverse(in_edge))
            }
            if (v0_edges.isEmpty()) continue

            // Sort the edges in clockwise order around "v0".
            val min_endpoint: VertexId = v0_edges.first().endpoint

            val v0VertexEdge = v0_edges.removeFirst()
            v0_edges.sortWith { a: VertexEdge, b: VertexEdge ->
                when {
                    a.endpoint == b.endpoint -> if(a.rank < b.rank) -1 else 1
                    a.endpoint == min_endpoint -> -1
                    b.endpoint == min_endpoint -> 1
                    !S2Predicates.orderedCCW(vertex(a.endpoint), vertex(b.endpoint), vertex(min_endpoint), vertex(v0)) -> -1
                    else -> 1
                }
            }
            v0_edges.add(0, v0VertexEdge)

            // Match incoming with outgoing edges.  We do this by keeping a stack of
            // unmatched incoming edges.  We also keep a stack of outgoing edges with
            // no previous incoming edge, and match these at the end by wrapping
            // around circularly to the start of the edge ordering.
            for (e in v0_edges) {
                when {
                    e.incoming -> e_in.add(in_edge_ids[e.index])
                    e_in.isNotEmpty() -> {
                        left_turn_map[e_in.last()] = e.index
                        e_in.removeLast()
                    }
                    else ->  e_out.add(e.index)  // Matched below.
                }
            }
            // Pair up additional edges using the fact that the ordering is circular.
            e_out.reverse()
            while (e_out.isNotEmpty() && e_in.isNotEmpty()) {
                left_turn_map[e_in.last()] = e_out.last()
                e_out.removeLast()
                e_in.removeLast()
            }
            // We only need to process unmatched incoming edges, since we are only
            // responsible for creating left turn map entries for those edges.
            if (e_in.isNotEmpty() && error.isOk()) {
                error.code = S2Error.BUILDER_EDGES_DO_NOT_FORM_LOOPS
                error.text = "Given edges do not form loops (indegree != outdegree)"
            }
            e_in.clear()
            e_out.clear()
            v0_edges.clear()
        }
        return error.isOk()
    }

    companion object {

        // Defines a value larger than any valid InputEdgeId.
        val kMaxInputEdgeId: InputEdgeId = Integer.MAX_VALUE

        // The following value of InputEdgeId means that an edge does not
        // corresponds to any input edge.
        val kNoInputEdgeId: InputEdgeId = kMaxInputEdgeId - 1

        // Given an edge (src, dst), returns the reverse edge (dst, src).
        fun reverse(e: Edge): Edge = Edge(e.second, e.first)

        // Rotates the edges of "loop" if necessary so that the edge(s) with the
        // largest input edge ids are last.  This ensures that when an output loop
        // is equivalent to an input loop, their cyclic edge orders are the same.
        // "min_input_ids" is the output of GetMinInputEdgeIds().
        fun canonicalizeLoopOrder(min_input_ids: List<InputEdgeId>, loop: MutableList<EdgeId>): Unit = TODO()

        // Sorts the given edge chains (i.e., loops or polylines) by the minimum
        // input edge id of each chains's first edge.  This ensures that when the
        // output consists of multiple loops or polylines, they are sorted in the
        // same order as they were provided in the input.
        fun canonicalizeVectorOrder(min_input_ids: List<InputEdgeId>, chains: List<List<EdgeId>>): Unit = TODO()

        ////////////////////////////////////////////////////////////////////////
        //////////////// Helper Functions for Creating Graphs //////////////////

        // Given an unsorted collection of edges, transform them according to the
        // given set of GraphOptions.  This includes actions such as discarding
        // degenerate edges; merging duplicate edges; and canonicalizing sibling
        // edge pairs in several possible ways (e.g. discarding or creating them).
        // The output is suitable for passing to the Graph constructor.
        //
        // If options.edge_type() == EdgeType::UNDIRECTED, then all input edges
        // should already have been transformed into a pair of directed edges.
        //
        // "input_ids" is a vector of the same length as "edges" that indicates
        // which input edges were snapped to each edge.  This vector is also updated
        // appropriately as edges are discarded, merged, etc.
        //
        // Note that "options" may be modified by this method: in particular, the
        // edge_type() can be changed if sibling_pairs() is CREATE or REQUIRE (see
        // the description of S2Builder::GraphOptions).
        fun processEdges(options: GraphOptions, edges: List<Edge>, input_ids: List<InputEdgeIdSetId>, id_set_lexicon: IdSetLexicon, error: S2Error): Unit = TODO()

        // Given a set of vertices and edges, removes all vertices that do not have
        // any edges and returned the new, minimal set of vertices.  Also updates
        // each edge in "edges" to correspond to the new vertex numbering.  (Note
        // that this method does *not* merge duplicate vertices, it simply removes
        // vertices of degree zero.)
        //
        // The new vertex ordering is a subsequence of the original ordering,
        // therefore if the edges were lexicographically sorted before calling this
        // method then they will still be sorted after calling this method.
        //
        // The extra argument "tmp" points to temporary storage used by this method.
        // All calls to this method from a single thread can reuse the same
        // temporary storage.  It should initially point to an empty vector.  This
        // can make a big difference to efficiency when this method is called many
        // times (e.g. to extract the vertices for different layers), since the
        // incremental running time for each layer becomes O(edges.size()) rather
        // than O(vertices.size() + edges.size()).
        fun FilterVertices(vertices: List<S2Point>, edges: List<Edge>, tmp: List<VertexId>): List<S2Point> = TODO()

        // A comparison function that allows stable sorting with std::sort (which is
        // fast but not stable).  It breaks ties between equal edges by comparing
        // their edge ids.
        fun stableLessThan(a: Edge, b: Edge, ai: EdgeId, bi: EdgeId): Boolean {
            // The following is simpler but the compiler (2016) doesn't optimize it as
            // well as it should:
            //   return make_pair(a, ai) < make_pair(b, bi);
            if (a.first < b.first) return true
            if (b.first < a.first) return false
            if (a.second < b.second) return true
            if (b.second < a.second) return false
            return ai < bi  // Stable sort.
        }

        // Given a set of duplicate outgoing edges (v0, v1) and a set of duplicate
        // incoming edges (v1, v0), this method assigns each edge an integer "rank" so
        // that the edges are sorted in a consistent order with respect to their
        // orderings around "v0" and "v1".  Usually there is just one edge, in which
        // case this is easy.  Sometimes there is one edge in each direction, in which
        // case the outgoing edge is always ordered before the incoming edge.
        //
        // In general, we allow any number of duplicate edges in each direction, in
        // which case outgoing edges are interleaved with incoming edges so as to
        // create as many degenerate (two-edge) loops as possible.  In order to get a
        // consistent ordering around "v0" and "v1", we move forwards through the list
        // of outgoing edges and backwards through the list of incoming edges.  If
        // there are more incoming edges, they go at the beginning of the ordering,
        // while if there are more outgoing edges then they go at the end.
        //
        // For example, suppose there are 2 edges "a,b" from "v0" to "v1", and 4 edges
        // "w,x,y,z" from "v1" to "v0".  Using lower/upper case letters to represent
        // incoming/outgoing edges, the clockwise ordering around v0 would be zyAxBw,
        // and the clockwise ordering around v1 would be WbXaYZ.  (Try making a
        // diagram with each edge as a separate arc.)
        private fun addVertexEdges(out_begin: EdgeId, out_end: EdgeId, in_begin: EdgeId, in_end: EdgeId, v1: VertexId, v0_edges: MutableList<VertexEdge>) {
            var rank = 0
            var in_end = in_end
            var out_begin = out_begin
            // Any extra incoming edges go at the beginning of the ordering.
            while (in_end - in_begin > out_end - out_begin) {
                v0_edges.add(VertexEdge(true, --in_end, v1, rank++))
            }
            // Next we interleave as many outgoing and incoming edges as possible.
            while (in_end > in_begin) {
                v0_edges.add(VertexEdge(false, out_begin++, v1, rank++))
                v0_edges.add(VertexEdge(true, --in_end, v1, rank++))
            }
            // Any extra outgoing edges to at the end of the ordering.
            while (out_end > out_begin) {
                v0_edges.add(VertexEdge(false, out_begin++, v1, rank++))
            }
        }
    }

    // A helper class for VertexOutMap that represents the outgoing edges
    // from a given vertex.
    data class VertexOutEdges(val beginIdx: Int, val endIdx: Int) {

        fun size(): Int {
            return endIdx - beginIdx
        }

    }
/*
    // A helper class for VertexOutMap that represents the outgoing edge *ids*
    // from a given vertex.
    class VertexOutEdgeIds
        : public std::iterator<std::forward_iterator_tag, EdgeId>
    {
        public:
        // An iterator over a range of edge ids (like boost::counting_iterator).
        class Iterator {
            public :
            explicit Iterator(EdgeId id) : id_(id)
            {}
            const EdgeId& operator *()
            const { return id_; }
            Iterator& operator ++()
            { ++id_; return * this; }
            Iterator operator ++(int)
            { return Iterator(id_++); }
            size_t operator -(const Iterator& x)
            const { return id_ - x.id_; }
            bool operator ==(const Iterator& x)
            const { return id_ == x.id_; }
            bool operator !=(const Iterator& x)
            const { return id_ != x.id_; }

            private :
            EdgeId id_;
        };
        Iterator begin () const { return Iterator(begin_); }
        Iterator end () const { return Iterator(end_); }
        size_t size () const { return end_ - begin_; }

        private:
        friend class VertexOutMap;
        VertexOutEdgeIds(EdgeId begin, EdgeId end);
        EdgeId begin_, end_;
    };



    // A class that maps vertices to their outgoing edge ids.  Example usage:
    //   VertexOutMap out(g);
    //   for (Graph::EdgeId e : out.edge_ids(v)) { ... }
    //   for (const Graph::Edge& edge : out.edges(v)) { ... }
    inner class VertexOutMap() {

        private val edgeBegins: List<EdgeId>

        // Return the edges (or edge ids) between a specific pair of vertices.
        fun edges(v: VertexId): VertexOutEdges {
            return VertexOutEdges(edgeBegins[v], edgeBegins[v + 1])
        }

        fun edges(v0: VertexId, v1: VertexId): VertexOutEdges {
            TODO()
            // auto range = std ::equal_range(edges_->data()+edge_begins_[v0],
            //        edges_->data()+edge_begins_[v0+1],
            //        Edge(v0, v1));
            //        return VertexOutEdges(range.first, range.second);
        }

        fun edge_ids(v0: VertexId, v1: VertexId): VertexOutEdgeIds {
            TODO()
//            auto range = std ::equal_range(edges_->data()+edge_begins_[v0],
//            edges_->data()+edge_begins_[v0+1],
//            Edge(v0, v1));
//            return VertexOutEdgeIds(
//                    static_cast < S2Builder::Graph::EdgeId > (range.first - edges_->data()),
//            static_cast < S2Builder::Graph::EdgeId > (range.second - edges_->data()));
        }

        fun degree(v: VertexId): Int = edge_ids(v).size()

        fun edge_ids(v: VertexId): VertexOutEdgeIds = VertexOutEdgeIds(edgeBegins[v], edgeBegins[v + 1])


        private :
        const std::vector<Edge>* edges_;
        std::vector<EdgeId> edge_begins_;
    }

    // A helper class for VertexInMap that represents the incoming edge *ids*
    // to a given vertex.
    class VertexInEdgeIds {
        public :
        const EdgeId* begin()
        const { return begin_; }
        const EdgeId* end()
        const { return end_; }
        size_t size()
        const { return end_ - begin_; }

        private :
        friend
        class VertexInMap;
        VertexInEdgeIds(const EdgeId* begin, const EdgeId* end);
        const EdgeId* begin_;
        const EdgeId* end_;
    };

    // A class that maps vertices to their incoming edge ids.  Example usage:
    //   VertexInMap in(g);
    //   for (Graph::EdgeId e : in.edge_ids(v)) { ... }
    class VertexInMap {
        public :
        VertexInMap() = default;
        explicit VertexInMap(const Graph& g)
        { Init(g); }
        void Init(const Graph& g);

        int degree(VertexId v) const ;
        VertexInEdgeIds edge_ids(VertexId v) const ;

        // Returns a sorted vector of all incoming edges (see GetInEdgeIds).
        const std::vector<EdgeId>& in_edge_ids()
        const { return in_edge_ids_; }

        private :
        std::vector<EdgeId> in_edge_ids_;
        std::vector<EdgeId> in_edge_begins_;
        VertexInMap(const VertexInMap&) = delete;
        void operator =(const VertexInMap&) = delete;
    }

    // Convenience class to return the set of labels associated with a given
    // graph edge.  Note that due to snapping, one graph edge may correspond to
    // several different input edges and will have all of their labels.
    // This class is the preferred way to retrieve edge labels.
    //
    // The reason this is a class rather than a graph method is because for
    // undirected edges, we need to fetch the labels associated with both
    // siblings.  This is because only the original edge of the sibling pair has
    // labels; the automatically generated sibling edge does not.
    class LabelFetcher {
        public :
        LabelFetcher() = default;
        LabelFetcher(const Graph& g, EdgeType edge_type)
        { Init(g, edge_type); }

        // Prepares to fetch labels associated with the given edge type.  For
        // EdgeType::UNDIRECTED, labels associated with both edges of the sibling
        // pair will be returned.  "edge_type" is a parameter (rather than using
        // g.options().edge_type()) so that clients can explicitly control whether
        // labels from one or both siblings are returned.
        void Init(const Graph& g, EdgeType edge_type);

        // Returns the set of labels associated with edge "e" (and also the labels
        // associated with the sibling of "e" if edge_type() is UNDIRECTED).
        // Labels are sorted and duplicate labels are automatically removed.
        //
        // This method uses an output parameter rather than returning by value in
        // order to avoid allocating a new vector on every call to this method.
        void Fetch(EdgeId e, std::vector<S2Builder::Label>* labels);

        private :
        const Graph* g_;
        EdgeType edge_type_;
        std::vector<EdgeId> sibling_map_;
    };
 */
    // Indicates whether loops should be simple cycles (no repeated vertices) or
    // circuits (which allow repeated vertices but not repeated edges).  In
    // terms of how the loops are built, this corresponds to closing off a loop
    // at the first repeated vertex vs. the first repeated edge.
    enum class LoopType { SIMPLE, CIRCUIT };

    // Builds loops from a set of directed edges, turning left at each vertex
    // until a repeated edge is found (i.e., LoopType::CIRCUIT).  The loops are
    // further grouped into connected components, where each component consists
    // of one or more loops connected by shared vertices.
    //
    // This method is used to build polygon meshes from directed or undirected
    // input edges.  To convert the output of this method into a mesh, the
    // client must determine how the loops in different components are related
    // to each other: for example, several loops from different components may
    // bound the same region on the sphere, in which case all of those loops are
    // combined into a single polygon.  (See s2shapeutil::BuildPolygonBoundaries
    // and s2builderutil::LaxPolygonVectorLayer for details.)
    //
    // Note that loops may include both edges of a sibling pair.  When several
    // such edges are connected in a chain or a spanning tree, they form a
    // zero-area "filament".  The entire loop may be a filament (i.e., a
    // degenerate loop with an empty interior), or the loop may have have
    // non-empty interior with several filaments that extend inside it, or the
    // loop may consist of several "holes" connected by filaments.  These
    // filaments do not change the interior of any loop, so if you are only
    // interested in point containment then they can safely be removed by
    // setting the "degenerate_boundaries" parameter to DISCARD.  (They can't be
    // removed by setting (options.sibling_pairs() == DISCARD) because the two
    // siblings might belong to different polygons of the mesh.)  Note that you
    // can prevent multiple copies of sibling pairs by specifying
    // options.duplicate_edges() == MERGE.
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
    //           (but requires DISCARD if degenerate_boundaries == DISCARD)
    // REQUIRES: options.sibling_pairs() == { REQUIRE, CREATE }
    //           [i.e., every edge must have a sibling edge]
    enum class DegenerateBoundaries { DISCARD, KEEP }

    fun getDirectedComponents(degenerate_boundaries: DegenerateBoundaries, components: MutableList<DirectedComponent>, error: S2Error): Boolean = TODO()

    // Builds loops from a set of undirected edges, turning left at each vertex
    // until either a repeated vertex (for LoopType::SIMPLE) or a repeated edge
    // (for LoopType::CIRCUIT) is found.  The loops are further grouped into
    // "components" such that all the loops in a component are connected by
    // shared vertices.  Finally, the loops in each component are divided into
    // two "complements" such that every edge in one complement is the sibling
    // of an edge in the other complement.  This corresponds to the fact that
    // given any set of non-crossing undirected loops, there are exactly two
    // possible interpretations of the region that those loops represent (where
    // one possibility is the complement of the other).  This method does not
    // attempt to resolve this ambiguity, but instead returns both possibilities
    // for each connected component and lets the client choose among them.
    //
    // This method is used to build single polygons.  (Use GetDirectedComponents
    // to build polygon meshes, even when the input edges are undirected.)  To
    // convert the output of this method into a polygon, the client must choose
    // one complement from each component such that the entire set of loops is
    // oriented consistently (i.e., they define a region such that the interior
    // of the region is always on the left).  The non-chosen complements form
    // another set of loops that are also oriented consistently but represent
    // the complementary region on the sphere.  Finally, the client needs to
    // choose one of these two sets of loops based on heuristics (e.g., the area
    // of each region), since both sets of loops are equally valid
    // interpretations of the input.
    //
    // Each loop is represented as a sequence of edges.  The edge ordering and
    // loop ordering are automatically canonicalized in order to preserve the
    // input ordering as much as possible.  Loops are non-crossing provided that
    // the graph contains no crossing edges.  If some edges cannot be turned
    // into loops, returns false and sets "error" appropriately.
    //
    // REQUIRES: options.degenerate_edges() == { DISCARD, DISCARD_EXCESS }
    // REQUIRES: options.edge_type() == UNDIRECTED
    // REQUIRES: options.siblings_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
    //           [since REQUIRE, CREATE convert the edge_type() to DIRECTED]
    fun getUndirectedComponents(loop_type: LoopType, components: MutableList<UndirectedComponent>, error: S2Error): Boolean = TODO()

    // Indicates whether polylines should be "paths" (which don't allow
    // duplicate vertices, except possibly the first and last vertex) or
    // "walks" (which allow duplicate vertices and edges).
    enum class PolylineType { PATH, WALK };

    // Builds polylines from a set of edges.  If "polyline_type" is PATH, then
    // only vertices of indegree and outdegree 1 (or degree 2 in the case of
    // undirected edges) will appear in the interior of polylines.  This
    // essentially generates one polyline for each edge chain in the graph.  If
    // "polyline_type" is WALK, then polylines may pass through the same vertex
    // or even the same edge multiple times (if duplicate edges are present),
    // and each polyline will be as long as possible.  This option is useful for
    // reconstructing a polyline that has been snapped to a lower resolution,
    // since snapping can cause edges to become identical.
    //
    // This method attempts to preserve the input edge ordering in order to
    // implement idempotency, even when there are repeated edges or loops.  This
    // is true whether directed or undirected edges are used.  Degenerate edges
    // are also handled appropriately.
    //
    // REQUIRES: options.sibling_pairs() == { DISCARD, DISCARD_EXCESS, KEEP }
    fun getPolylines(polyline_type: PolylineType): List<EdgePolyline> = TODO()

}
