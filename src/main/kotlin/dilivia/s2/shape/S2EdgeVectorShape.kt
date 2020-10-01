package dilivia.s2.shape

import dilivia.s2.Assertions
import dilivia.s2.S2Point


// S2EdgeVectorShape is an S2Shape representing an arbitrary set of edges.  It
// is mainly used for testing, but it can also be useful if you have, say, a
// collection of polylines and don't care about memory efficiency (since this
// class would store most of the vertices twice).
//
// Note that if you already have data stored in an S2Loop, S2Polyline, or
// S2Polygon, then you would be better off using the "Shape" class defined
// within those classes (e.g., S2Loop::Shape).  Similarly, if the vertex data
// is stored in your own data structures, you can easily write your own
// subclass of S2Shape that points to the existing vertex data rather than
// copying it.
class S2EdgeVectorShape(id: Int = -1, edges: List<Edge> = emptyList()) : S2Shape(id) {

    private val edges = edges.toMutableList()

  // Creates an S2EdgeVectorShape containing a single edge.
  constructor(a: S2Point, b: S2Point): this(edges = listOf(Edge(a, b)))

  // Adds an edge to the vector.
  //
  // IMPORTANT: This method should only be called *before* adding the
  // S2EdgeVectorShape to an S2ShapeIndex.  S2Shapes can only be modified by
  // removing them from the index, making changes, and adding them back again.
  fun add(a: S2Point, b: S2Point) {
    edges.add(Edge(a, b))
  }

  // S2Shape interface:

    override val numEdges: Int
        get() = edges.size

    override fun edge(edgeId: Int): Edge = edges[edgeId]

    override val dimension: Int = 1

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int
        get() = edges.size

    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        Assertions.assertEQ(offset, 0)
        return edge(chainId)
    }

    override fun chainPosition(edgeId: Int): ChainPosition  = ChainPosition(edgeId, 0)

    override val typeTag: TypeTag = kNoTypeTag

}
