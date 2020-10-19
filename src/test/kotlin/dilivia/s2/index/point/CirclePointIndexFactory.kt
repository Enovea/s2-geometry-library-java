package dilivia.s2.index.point

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.region.S2Cap

// Generates points that are regularly spaced along a circle.  Points along a
// circle are nearly the worst case for distance calculations, since many
// points are nearly equidistant from any query point that is not immediately
// adjacent to the circle.
class CirclePointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val points = S2GeometryTestCase.makeRegularPoints(indexCap.center, indexCap.radius(), numPoints)
        points.forEachIndexed { i, point -> index.add(point, i) }
    }
}
