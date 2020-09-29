package dilivia.s2.shape

import dilivia.s2.S1Angle
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.region.S2Loop
import dilivia.s2.shape.S2ShapeUtil.containsBruteForce

/*
class S2ShapeUtilTest : S2GeometryTestCase() {

    fun testContainsBruteForceNoInterior() {
        // Defines a polyline that almost entirely encloses the point 0:0.
        val polyline = makeLaxPolyline("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        assertFalse(S2ShapeUtil.containsBruteForce(polyline, makePoint("0:0")))
    }

    fun testContainsBruteForceContainsReferencePoint() {
        // Checks that ContainsBruteForce agrees with GetReferencePoint.
        val polygon = makeLaxPolygon("0:0, 0:1, 1:-1, -1:-1, -1e9:1");
        val ref = polygon.getReferencePoint()
        assertEquals(ref.contained, containsBruteForce(polygon, ref.point))
    }

    fun testContainsBruteForceConsistentWithS2Loop() {
        // Checks that ContainsBruteForce agrees with S2Loop::Contains().
        val loop = S2Loop.makeRegularLoop(makePoint("89:-179"), S1Angle.degrees(10), 100)
        val shape = S2Loop.Shape(0 , loop)
        for (i in 0 until loop.numVertices()) {
            assertEquals(loop.contains(loop.vertex(i)), S2ShapeUtil.containsBruteForce(shape, loop.vertex(i)))
        }
    }


}

 */