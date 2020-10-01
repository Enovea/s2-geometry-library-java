package dilivia.s2.index

object S2CountEdges {

    // Returns the total number of edges in all indexed shapes.  This method takes
    // time linear in the number of shapes.
    fun countEdges(index: S2ShapeIndex): Int = countEdgesUpTo(index, Int.MAX_VALUE)

    // Like CountEdges(), but stops once "max_edges" edges have been found (in
    // which case the current running total is returned).
    fun countEdgesUpTo(index: S2ShapeIndex, max_edges: Int): Int {
        var num_edges = 0
        val shapeIter = index.begin()
        shapeIter.asSequence().filterNotNull().forEach { shape ->
            num_edges += shape.numEdges
            if (num_edges >= max_edges) return num_edges
        }
        return num_edges;
    }


}