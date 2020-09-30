package dilivia.s2.shape

// S2ClippedShape represents the part of a shape that intersects an S2Cell.
// It consists of the set of edge ids that intersect that cell, and a boolean
// indicating whether the center of the cell is inside the shape (for shapes
// that have an interior).
//
// Note that the edges themselves are not clipped; we always use the original
// edges for intersection tests so that the results will be the same as the
// original shape.
// @property shapeId The shape id of the clipped shape.
data class S2ClippedShape(val shapeId: Int, val edges: List<Int>, val containsCenter: Boolean) {

    // Returns true if the center of the S2CellId is inside the shape.  Returns
    // false for shapes that do not have an interior.
    fun containsCenter(): Boolean = containsCenter

    // The number of edges that intersect the S2CellId.
    fun numEdges(): Int = edges.size

    // Returns the edge id of the given edge in this clipped shape.  Edges are
    // sorted in increasing order of edge id.
    //
    // REQUIRES: 0 <= i < num_edges()
    fun edge(i: Int): Int = edges[i]

    // Returns true if the clipped shape contains the given edge id.
    fun containsEdge(id: Int): Boolean {
        // Linear search is fast because the number of edges per shape is typically
        // very small (less than 10).
        for (e in 0 until numEdges()) {
            if (edge(e) == id) return true
        }
        return false
    }

}