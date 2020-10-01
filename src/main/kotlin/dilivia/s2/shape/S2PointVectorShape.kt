package dilivia.s2.shape

import dilivia.s2.Assertions
import dilivia.s2.S2Point

// S2PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is reprsented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.
class S2PointVectorShape(id: Int = -1, val points: List<S2Point> = emptyList()) : S2Shape(id) {

  fun numPoints(): Int = points.size

  fun point(i: Int): S2Point = points[i]

  // S2Shape interface:

    override val numEdges: Int
        get() = numPoints()

    override fun edge(edgeId: Int): Edge = Edge(points[edgeId], points[edgeId])

    override val dimension: Int = 0

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int
        get() = numPoints()

    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)

    override fun chainEdge(chainId: Int, offset: Int): Edge {
        Assertions.assertEQ(offset, 0)
        return Edge(points[chainId], points[chainId])
    }

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    override val typeTag: TypeTag = TypeTags.kPointVectorTypeTag

};
