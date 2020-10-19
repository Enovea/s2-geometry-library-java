package dilivia.s2.index.point

import dilivia.s2.S1ChordAngle
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser
import dilivia.s2.index.S2MaxDistance
import dilivia.s2.index.S2MaxDistanceFactory
import dilivia.s2.index.S2MaxDistancePointTarget

//
// This file contains some basic tests of the templating support.  Testing of
// the actual algorithms is in s2closest_point_query_test.cc.


// This is a proof-of-concept prototype of a possible S2FurthestPointQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
typealias FurthestPointQuery = S2ClosestPointQueryBase<S2MaxDistance, Int>

class FurthestPointTarget(point: S2Point) : S2MaxDistancePointTarget(point) {

    override fun maxBruteForceIndexSize(): Int = 10

}

class S2ClosestPointQueryBaseTest : S2GeometryTestCase() {

    fun testS2ClosestPointQueryBaseMaxDistance() {
        val index = S2PointIndex<Int>()
      val points = S2TextParser.parsePoints("0:0, 1:0, 2:0, 3:0")
        for (i in points.indices) {
            index.add(points[i], i);
        }
        val query = FurthestPointQuery(S2MaxDistanceFactory, index)
        val options = S2ClosestPointQueryBase.Options(S2MaxDistanceFactory)
        options.setMaxResult(1)
        val target = FurthestPointTarget(S2TextParser.makePoint("4:0"))
        val results = query . findClosestPoints (target, options)
        assertEquals(1, results.size);
        assertEquals(points[0], results[0].point());
        assertEquals(0, results[0].data());
        assertDoubleNear(4.0, S1ChordAngle(results[0].distance.value).toAngle().degrees(), 1e-13)
    }

}
