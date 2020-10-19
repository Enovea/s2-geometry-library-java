package dilivia.s2.index.point

import dilivia.s2.region.S2Cap

// An abstract class that adds points to an S2PointIndex for benchmarking.
interface PointIndexFactory {

    // Requests that approximately "num_points" points located within the given
    // S2Cap bound should be added to "index".
    fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>)

}
