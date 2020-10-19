package dilivia.s2.index.shape

import dilivia.s2.region.S2Cap

// An abstract class that adds edges to a MutableS2ShapeIndex for benchmarking.
interface ShapeIndexFactory {

    // Requests that approximately "num_edges" edges located within the given
    // S2Cap bound should be added to "index".
    fun addEdges(index_cap: S2Cap, num_edges: Int, index: MutableS2ShapeIndex)

}
