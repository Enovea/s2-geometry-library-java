package dilivia.s2.shape

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32, int32> except that it has named fields.
// It should be passed and returned by value.
data class ShapeEdgeId(val shapeId: Int = -1, val edgeId: Int = -1): Comparable<ShapeEdgeId> {

    override fun compareTo(other: ShapeEdgeId): Int {
        val shapeIdComparison = shapeId.compareTo(other.shapeId)
        return if (shapeIdComparison != 0) shapeIdComparison else edgeId.compareTo(other.edgeId)
    }

}
