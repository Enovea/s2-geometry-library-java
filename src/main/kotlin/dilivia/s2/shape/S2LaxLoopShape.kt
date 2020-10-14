package dilivia.s2.shape

import dilivia.s2.Assertions
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.Assertions.assertLT
import dilivia.s2.S2Point
import dilivia.s2.region.S2Loop
import kotlin.math.min

//
// This file defines various S2Shape types representing loops:
//
// S2LaxLoopShape
//   - like S2Loop::Shape but allows duplicate vertices & edges, more compact
//     representation, and faster to initialize.
//
// S2LaxClosedPolylineShape
//   - like S2LaxLoopShape, but defines a loop that does not have an interior
//     (a closed polyline).
//
// S2VertexIdLaxLoopShape
//   - like S2LaxLoopShape, but vertices are specified as indices into an
//     existing vertex array.


// S2LaxLoopShape represents a closed loop of edges surrounding an interior
// region.  It is similar to S2Loop::Shape except that this class allows
// duplicate vertices and edges.  Loops may have any number of vertices,
// including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
// consisting of a single point.)
//
// Note that S2LaxLoopShape is faster to initialize and more compact than
// S2Loop::Shape, but does not support the same operations as S2Loop.
open class S2LaxLoopShape : S2Shape {

    // For clients that have many small loops, we save some memory by
    // representing the vertices as an array rather than using std::vector.
    private lateinit var vertices: Array<S2Point>

    // Constructs an empty loop.
    constructor() {
        vertices = emptyArray()
    }

    // Constructs an S2LaxLoopShape with the given vertices.
    constructor(vertices: List<S2Point>) {
        init(vertices)
    }

    // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
    constructor(loop: S2Loop) {
        init(loop)
    }

    // Initializes an S2LaxLoopShape with the given vertices.
    fun init(vertices: List<S2Point>) {
        this.vertices = vertices.toTypedArray()
    }

    // Initializes an S2LaxLoopShape from the given S2Loop, by copying its data.
    //
    // REQUIRES: !loop->is_full()
    //           [Use S2LaxPolygonShape if you need to represent a full loop.]
    fun init(loop: S2Loop) {
        Assertions.assert({ !loop.isFull() }, { "Full loops not supported; use S2LaxPolygonShape" })
        if (loop.isEmpty()) {
            vertices = emptyArray()
        } else {
            vertices = loop.verticesSpan().points.toTypedArray()
        }
    }

    fun numVertices(): Int = vertices.size
    fun vertex(i: Int): S2Point = vertices[i]

    override val numEdges: Int
        get() = numVertices()

    override fun edge(edgeId: Int): Edge {
        Assertions.assertLT(edgeId, numEdges)
        var e = edgeId + 1
        if (e == numVertices()) e = 0
        return Edge(vertices[edgeId], vertices[e])
    }

    override val dimension: Int = 2

    override fun getReferencePoint(): ReferencePoint = S2Shape.getReferencePoint(this)

    override val numChains: Int
        get() = min(1, numVertices())

    override fun chain(chain_id: Int): Chain = Chain(0, numVertices())

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        assertEQ(chainId, 0)
        assertLT(offset, numEdges)
        val k = if(offset + 1 == numVertices()) 0 else offset + 1
        return Edge(vertices[offset], vertices[k])
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

}

// S2LaxClosedPolylineShape is like S2LaxPolylineShape except that the last
// vertex is implicitly joined to the first.  It is also like S2LaxLoopShape
// except that it does not have an interior (which makes it more efficient to
// index).
class S2LaxClosedPolylineShape : S2LaxLoopShape {

    // Constructs an empty loop.
    constructor(): super()

    // Constructs an S2LaxLoopShape with the given vertices.
    constructor(vertices: List<S2Point>): super(vertices)

    // Constructs an S2LaxLoopShape from the given S2Loop, by copying its data.
    constructor(loop: S2Loop): super(loop)

    override val dimension: Int = 1

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

}

/*
// S2VertexIdLaxLoopShape is just like S2LaxLoopShape, except that vertices are
// specified as indices into a vertex array.  This representation can be more
// compact when many loops are arranged in a mesh structure.
class S2VertexIdLaxLoopShape : S2Shape() {
    public:
    // Constructs an empty loop.
    S2VertexIdLaxLoopShape() : num_vertices_(0) {}

    // Constructs the shape from the given vertex array and indices.
    // "vertex_ids" is a vector of indices into "vertex_array".
    //
    // ENSURES:  loop->vertex(i) == (*vertex_array)[vertex_ids[i]]
    // REQUIRES: "vertex_array" persists for the lifetime of this object.
    explicit S2VertexIdLaxLoopShape (const std ::vector<int32>& vertex_ids,
    const S2Point * vertex_array);

    // Initializes the shape from the given vertex array and indices.
    // "vertex_ids" is a vector of indices into "vertex_array".
    void Init (const std ::vector<int32>& vertex_ids,
    const S2Point * vertex_array);

    // Returns the number of vertices in the loop.
    int num_vertices () const { return num_vertices_; }
    int32 vertex_id (int i) const { return vertex_ids_[i]; }
    const S2Point & vertex (int i) const { return vertex_array_[vertex_id(i)]; }

    // S2Shape interface:
    int num_edges () const final { return num_vertices(); }
    Edge edge (int e) const final;
    int dimension () const final { return 2; }
    ReferencePoint GetReferencePoint () const final;
    int num_chains () const final { return std::min(1, num_vertices_); }
    Chain chain (int i) const final { return Chain(0, num_vertices_); }
    Edge chain_edge (int i, int j) const final;
    ChainPosition chain_position (int e) const final {
        return ChainPosition(0, e);
    }

    private:
    int32 num_vertices_;
    std::unique_ptr < int32[] > vertex_ids_;
    const S2Point * vertex_array_;
};
*/
