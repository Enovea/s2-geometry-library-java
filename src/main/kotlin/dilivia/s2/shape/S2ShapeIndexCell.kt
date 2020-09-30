package dilivia.s2.shape

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
class S2ShapeIndexCell() {

    private val shapes = mutableListOf<S2ClippedShape>()

    // Returns the number of clipped shapes in this cell.
    fun numClipped(): Int {
        return shapes.size
    }

    // Returns the clipped shape at the given index.  Shapes are kept sorted in
    // increasing order of shape id.
    //
    // REQUIRES: 0 <= i < num_clipped()
    fun clipped(i: Int): S2ClippedShape {
        return shapes[i]
    }

    // Returns a pointer to the clipped shape corresponding to the given shape,
    // or nullptr if the shape does not intersect this cell.
    fun findClipped(shape: S2Shape): S2ClippedShape? = findClipped(shape.id)
    fun findClipped(shape_id: Int): S2ClippedShape? {
        // Linear search is fine because the number of shapes per cell is typically
        // very small (most often 1), and is large only for pathological inputs
        // (e.g. very deeply nested loops).
        for (s in shapes) {
            if (s.shapeId == shape_id) return s
        }
        return null
    }

    // Convenience method that returns the total number of edges in all clipped
    // shapes.
    fun num_edges(): Int {
        var n = 0
        for (i in 0 until numClipped()) n += clipped(i).numEdges()
        return n
    }

    fun addClipped(shape: S2ClippedShape) {
        shapes.add(shape)
    }

    override fun toString(): String {
        return "S2ShapeIndexCell(shapes=$shapes)"
    }


}