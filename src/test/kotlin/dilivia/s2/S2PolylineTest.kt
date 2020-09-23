package dilivia.s2

import com.google.common.geometry.S2.M_PI
import mu.KotlinLogging
import kotlin.math.*

class S2PolylineTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger {}
    
    fun testBasic() {
        val vertices = emptyList<S2Point>()
        var empty = S2Polyline(vertices)
        assertEquals(S2LatLngRect.empty, empty.rectBound)
        empty = empty.reversed()
        assertEquals(0, empty.numVertices())

        val latlngs = listOf(
                S2LatLng.fromDegrees(0, 0),
                S2LatLng.fromDegrees(0, 90),
                S2LatLng.fromDegrees(0, 180)
        )
        var semiEquator = S2Polyline.fromLatLng(latlngs)
        assertTrue(S2Point.approxEquals(semiEquator.interpolate(0.5), S2Point(0, 1, 0)))
        semiEquator = semiEquator.reversed()
        assertEquals(S2Point(1, 0, 0), semiEquator.vertex(2))
    }

    fun testGetLengthAndCentroid() {
        // Construct random great circles and divide them randomly into segments.
        // Then make sure that the length and centroid are correct.  Note that
        // because of the way the centroid is computed, it does not matter how
        // we split the great circle into segments.

        repeat(100) {
            // Choose a coordinate frame for the great circle.
            val (x, y, _) = S2Random.randomFrame()

            val vertices = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                val p = cos(theta) * x + sin(theta) * y
                if (vertices.isEmpty() || p != vertices.last())
                    vertices.add(p)
                theta += S2Random.randomDouble().pow(10.0)
            }
            // Close the circle.
            vertices.add(vertices[0])
            val line = S2Polyline(vertices)
            val length = line.getLength()
            assertTrue(abs(length.radians - 2 * M_PI) <= 2e-14)
            val centroid = line.getCentroid()
            assertTrue(centroid.norm() <= 2e-14)
        }
    }

    fun testMayIntersect() {
        val vertices = listOf(
                S2Point(1.0, -1.1, 0.8).normalize(),
                S2Point(1.0, -0.8, 1.1).normalize()
        )
        val line = S2Polyline(vertices)
        for (face in 0..5) {
            val cell = S2Cell.fromFace(face)
            assertEquals((face and 1) == 0, line.mayIntersect(cell))
        }
    }

    fun testInterpolate() {
        var vertices = listOf(
                S2Point(1, 0, 0),
                S2Point(0, 1, 0),
                S2Point(0, 1, 1).normalize(),
                S2Point(0, 0, 1)
        )
        val line = S2Polyline(vertices)
        assertEquals(vertices[0], line.interpolate(-0.1))
        assertTrue(S2Point.approxEquals(line.interpolate(0.1), S2Point(1.0, tan(0.2 * M_PI / 2), 0.0).normalize()))
        assertTrue(S2Point.approxEquals(line.interpolate(0.25), S2Point(1, 1, 0).normalize()))
        assertEquals(vertices[1], line.interpolate(0.5))
        assertTrue(S2Point.approxEquals(vertices[2], line.interpolate(0.75)))
        val (vertex0, next_vertex0) = line.getSuffix(-0.1)
        assertEquals(vertices[0], vertex0)
        assertEquals(1, next_vertex0)
        val (vertex2, next_vertex2) = line.getSuffix(0.75)
        assertTrue(S2Point.approxEquals(vertices[2], vertex2))
        assertEquals(3, next_vertex2)
        val (vertex3, next_vertex3) = line.getSuffix(1.1)
        assertEquals(vertices[3], vertex3)
        assertEquals(4, next_vertex3)

        // Check the case where the interpolation fraction is so close to 1 that
        // the interpolated point is identical to the last vertex.
        vertices = listOf(
            S2Point(1, 1, 1).normalize(),
            S2Point(1.0, 1.0, 1 + 1e-15).normalize(),
            S2Point(1.0, 1.0, 1 + 2e-15).normalize()
        )
        val shortLine = S2Polyline(vertices)
        val (vertex, next_vertex) = shortLine.getSuffix(1.0 - 2e-16)
        assertEquals(vertices[2], vertex)
        assertEquals(3, next_vertex)
    }

    fun testUnInterpolate() {
        val vertices = mutableListOf( S2Point(1, 0, 0) )
        val pointLine = S2Polyline(vertices)
        assertEquals(0.0, pointLine.unInterpolate(S2Point(0, 1, 0), 1))

        vertices.add(S2Point(0, 1, 0))
        vertices.add(S2Point(0, 1, 1).normalize())
        vertices.add(S2Point(0, 0, 1))
        val line = S2Polyline(vertices)

        var interpolated: Pair<S2Point, Int>
        interpolated = line.getSuffix(-0.1)
        assertEquals(0.0, line.unInterpolate(interpolated.first, interpolated.second))
        interpolated = line.getSuffix(0.0)
        assertEquals(0.0, line.unInterpolate(interpolated.first, interpolated.second))
        interpolated = line.getSuffix(0.5)
        assertEquals(0.5, line.unInterpolate(interpolated.first, interpolated.second))
        interpolated = line.getSuffix(0.75)
        assertEquals(0.75, line.unInterpolate(interpolated.first, interpolated.second))
        interpolated = line.getSuffix(1.1)
        assertEquals(1.0, line.unInterpolate(interpolated.first, interpolated.second))

        // Check that the return value is clamped to 1.0.
        assertEquals(1.0, line.unInterpolate(S2Point(0, 1, 0), vertices.size))
    }

    fun testProject() {
        val latlngs = listOf(
            S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
            S2LatLng.fromDegrees(0, 2), S2LatLng.fromDegrees(1, 2)
        )
        val line = S2Polyline.fromLatLng(latlngs)


        var projected = line.project(S2LatLng.fromDegrees(0.5, -0.5).toPoint())
        assertTrue(S2Point.approxEquals(projected.first, S2LatLng.fromDegrees(0, 0).toPoint()))
        assertEquals(1, projected.second)
        projected = line.project(S2LatLng.fromDegrees(0.5, 0.5).toPoint())
        assertTrue(S2Point.approxEquals(projected.first, S2LatLng.fromDegrees(0.0, 0.5).toPoint()))
        assertEquals(1, projected.second)
        projected = line.project(S2LatLng.fromDegrees(0.5, 1.0).toPoint())
        assertTrue(S2Point.approxEquals(projected.first, S2LatLng.fromDegrees(0, 1).toPoint()))
        assertEquals(2, projected.second)
        projected = line.project(S2LatLng.fromDegrees(-0.5, 2.5).toPoint())
        assertTrue(S2Point.approxEquals(projected.first, S2LatLng.fromDegrees(0, 2).toPoint()))
        assertEquals(3, projected.second)
        projected = line.project(S2LatLng.fromDegrees(2, 2).toPoint())
        assertTrue(S2Point.approxEquals(projected.first, S2LatLng.fromDegrees(1, 2).toPoint()))
        assertEquals(4, projected.second)
    }

    fun testIsOnRight() {
        var latlngs = listOf(
            S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
            S2LatLng.fromDegrees(0, 2), S2LatLng.fromDegrees(1, 2)
        )
        val line = S2Polyline.fromLatLng(latlngs)

        assertTrue(line.isOnRight(S2LatLng.fromDegrees(-0.5, 0.5).toPoint()))
        assertFalse(line.isOnRight(S2LatLng.fromDegrees(0.5, -0.5).toPoint()))
        assertFalse(line.isOnRight(S2LatLng.fromDegrees(0.5, 0.5).toPoint()))
        assertFalse(line.isOnRight(S2LatLng.fromDegrees(0.5, 1.0).toPoint()))
        assertTrue(line.isOnRight(S2LatLng.fromDegrees(-0.5, 2.5).toPoint()))
        assertTrue(line.isOnRight(S2LatLng.fromDegrees(1.5, 2.5).toPoint()))

        // Explicitly test the case where the closest point is an interior vertex.
        latlngs = listOf(
            S2LatLng.fromDegrees(0, 0), S2LatLng.fromDegrees(0, 1),
            S2LatLng.fromDegrees(-1, 0)
        )
        val line2 = S2Polyline.fromLatLng(latlngs)

        // The points are chosen such that they are on different sides of the two
        // edges that the interior vertex is on.
        assertFalse(line2.isOnRight(S2LatLng.fromDegrees(-0.5, 5.0).toPoint()))
        assertFalse(line2.isOnRight(S2LatLng.fromDegrees(5.5, 5.0).toPoint()))
    }

    fun testIntersectsEmptyPolyline() {
        val line1 = makePolyline("1:1, 4:4")
        val emptyPolyline = S2Polyline()
        assertFalse(emptyPolyline.intersects(line1))
    }

    fun testIntersectsOnePointPolyline() {
        val line1 = makePolyline("1:1, 4:4")
        val line2 = makePolyline("1:1")
        assertFalse(line1.intersects(line2))
    }

    fun testIntersects() {
        val line1  = makePolyline("1:1, 4:4")
        val smallCrossing  = makePolyline("1:2, 2:1")
        val smallNoncrossing  = makePolyline("1:2, 2:3")
        val bigCrossing  = makePolyline("1:2, 2:3, 4:3")

        assertTrue(line1.intersects(smallCrossing))
        assertFalse(line1.intersects(smallNoncrossing))
        assertTrue(line1.intersects(bigCrossing))
    }

    fun testIntersectsAtVertex() {
        val line1  = makePolyline("1:1, 4:4, 4:6")
        val line2  = makePolyline("1:1, 1:2")
        val line3  = makePolyline("5:1, 4:4, 2:2")
        assertTrue(line1.intersects(line2))
        assertTrue(line1.intersects(line3))
    }

    fun testIntersectsVertexOnEdge() {
        val horizontalLeftToRight  = makePolyline("0:1, 0:3")
        val verticalBottomToTop  = makePolyline("-1:2, 0:2, 1:2")
        val horizontalRightToLeft  = makePolyline("0:3, 0:1")
        val verticalTopToBottom  = makePolyline("1:2, 0:2, -1:2")
        assertTrue(horizontalLeftToRight.intersects(verticalBottomToTop))
        assertTrue(horizontalLeftToRight.intersects(verticalTopToBottom))
        assertTrue(horizontalRightToLeft.intersects(verticalBottomToTop))
        assertTrue(horizontalRightToLeft.intersects(verticalTopToBottom))
    }

    fun checkSubsample(polyline_str: String, tolerance_degrees: Double, expected_str: String) {
        logger.info { """"$polyline_str, tolerance $tolerance_degrees""" }
        val polyline = makePolyline(polyline_str)
        val indices = polyline.subsampleVertices(S1Angle.degrees(tolerance_degrees))
        assertEquals(expected_str, indices.joinToString(","))
    }

    fun testSubsampleVerticesTrivialInputs() {
        // No vertices.
        checkSubsample("", 1.0, "")
        // One vertex.
        checkSubsample("0:1", 1.0, "0")
        // Two vertices.
        checkSubsample("10:10, 11:11", 5.0, "0,1")
        // Three points on a straight line.
        // In theory, zero tolerance should work, but in practice there are floating
        // point errors.
        checkSubsample("-1:0, 0:0, 1:0", 1e-15, "0,2")
        // Zero tolerance on a non-straight line.
        checkSubsample("-1:0, 0:0, 1:1", 0.0, "0,1,2")
        // Negative tolerance should return all vertices.
        checkSubsample("-1:0, 0:0, 1:1", -1.0, "0,1,2")
        // Non-zero tolerance with a straight line.
        checkSubsample("0:1, 0:2, 0:3, 0:4, 0:5", 1.0, "0,4")

        // And finally, verify that we still do something reasonable if the client
        // passes in an invalid polyline with two or more adjacent vertices.
        checkSubsample("0:1, 0:1, 0:1, 0:2", 0.0, "0,3")
    }

    fun testSubsampleVerticesSimpleExample() {
        val polyStr = "0:0, 0:1, -1:2, 0:3, 0:4, 1:4, 2:4.5, 3:4, 3.5:4, 4:4"
        checkSubsample(polyStr, 3.0, "0,9")
        checkSubsample(polyStr, 2.0, "0,6,9")
        checkSubsample(polyStr, 0.9, "0,2,6,9")
        checkSubsample(polyStr, 0.4, "0,1,2,3,4,6,9")
        checkSubsample(polyStr, 0.0, "0,1,2,3,4,5,6,7,8,9")
    }

    fun testSubsampleVerticesGuarantees() {
        // Check that duplicate vertices are never generated.
        checkSubsample("10:10, 12:12, 10:10", 5.0, "0")
        checkSubsample("0:0, 1:1, 0:0, 0:120, 0:130", 5.0, "0,3,4")

        // Check that points are not collapsed if they would create a line segment
        // longer than 90 degrees, and also that the code handles original polyline
        // segments longer than 90 degrees.
        checkSubsample("90:0, 50:180, 20:180, -20:180, -50:180, -90:0, 30:0, 90:0", 5.0, "0,2,4,5,6,7")

        // Check that the output polyline is parametrically equivalent and not just
        // geometrically equivalent, i.e. that backtracking is preserved.  The
        // algorithm achieves this by requiring that the points must be encountered
        // in increasing order of distance along each output segment, except for
        // points that are within "tolerance" of the first vertex of each segment.
        checkSubsample("10:10, 10:20, 10:30, 10:15, 10:40", 5.0, "0,2,3,4")
        checkSubsample("10:10, 10:20, 10:30, 10:10, 10:30, 10:40", 5.0, "0,2,3,5")
        checkSubsample("10:10, 12:12, 9:9, 10:20, 10:30", 5.0, "0,4")
    }


    fun testEquals(a_str: String, b_str: String, max_error: S1Angle): Boolean {
        val a = makePolyline(a_str)
        val b = makePolyline(b_str)
        return a.approxEquals(b, max_error)
    }

    fun testApproxEquals() {
        val degree = S1Angle.degrees(1)

        // Close lines, differences within max_error.
        assertTrue(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.5 * degree))

        // Close lines, differences outside max_error.
        assertFalse(testEquals("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.01 * degree))

        // Same line, but different number of vertices.
        assertFalse(testEquals("0:0, 0:10, 0:20", "0:0, 0:20", 0.1 * degree))

        // Same vertices, in different order.
        assertFalse(testEquals("0:0, 5:5, 0:10", "5:5, 0:10, 0:0", 0.1 * degree))
    }

    fun testS2PolylineShapeBasic() {
        val polyline  = makePolyline("0:0, 1:0, 1:1, 2:1")
        val shape = S2Polyline.Shape(polyline = polyline)
        assertEquals(polyline, shape.polyline)
        assertEquals(3, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(3, shape.chain(0).length)
        val edge2 = shape . edge (2)
        assertEquals(S2LatLng.fromDegrees(1, 1).toPoint(), edge2.v0)
        assertEquals(S2LatLng.fromDegrees(2, 1).toPoint(), edge2.v1)
        assertEquals(1, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testS2PolylineShapeEmptyPolyline() {
        val polyline = S2Polyline()
        val shape = S2Polyline.Shape(polyline = polyline)
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testNearlyCovers(a_str: String, b_str: String, max_error_degrees: Double, expect_b_covers_a: Boolean, expect_a_covers_b: Boolean) {
        fail("infinite loop")
        logger.info { "a=\"$a_str\", b=\"$b_str\", max error=$max_error_degrees" }
        val a = makePolyline(a_str)
        val b = makePolyline(b_str)
        val maxError = S1Angle.degrees(max_error_degrees)
        assertEquals(expect_b_covers_a, b.nearlyCovers(a, maxError))
        assertEquals(expect_a_covers_b, a.nearlyCovers(b, maxError))
    }

    fun testS2PolylineCoveringTestPolylineOverlapsSelf() {
        val pline = "1:1, 2:2, -1:10"
        testNearlyCovers(pline, pline, 1e-10, true, true)
    }

    fun testS2PolylineCoveringTestPolylineDoesNotOverlapReverse() {
        testNearlyCovers("1:1, 2:2, -1:10", "-1:10, 2:2, 1:1", 1e-10, false, false)
    }

    fun testS2PolylineCoveringTestPolylineOverlapsEquivalent() {
        // These two polylines trace the exact same polyline, but the second one uses
        // three points instead of two.
        testNearlyCovers("1:1, 2:1", "1:1, 1.5:1, 2:1", 1e-10, true, true)
    }

    fun testS2PolylineCoveringTestShortCoveredByLong() {
        // The second polyline is always within 0.001 degrees of the first polyline,
        // but the first polyline is too long to be covered by the second.
        testNearlyCovers(
                "-5:1, 10:1, 10:5, 5:10", "9:1, 9.9995:1, 10.0005:5", 1e-3, false, true)
    }

    fun testS2PolylineCoveringTestPartialOverlapOnly() {
        // These two polylines partially overlap each other, but neither fully
        // overlaps the other.
        testNearlyCovers("-5:1, 10:1", "0:1, 20:1", 1.0, false, false)
    }

    fun testS2PolylineCoveringTestShortBacktracking() {
        // Two lines that backtrack a bit (less than 1.5 degrees) on different edges.
        // A simple greedy matching algorithm would fail on this example.
        val t1 = "0:0, 0:2, 0:1, 0:4, 0:5"
        val t2 = "0:0, 0:2, 0:4, 0:3, 0:5"
        testNearlyCovers(t1, t2, 1.5, true, true)
        testNearlyCovers(t1, t2, 0.5, false, false)
    }

    fun testS2PolylineCoveringTestLongBacktracking() {
        // Two arcs with opposite direction do not overlap if the shorter arc is
        // longer than max_error, but do if the shorter arc is shorter than max-error.
        testNearlyCovers("5:1, -5:1", "1:1, 3:1", 1.0, false, false)
        testNearlyCovers("5:1, -5:1", "1:1, 3:1", 2.5, false, true)
    }

    fun testS2PolylineCoveringTestIsResilientToDuplicatePoints() {
        // S2Polyines are not generally supposed to contain adjacent, identical
        // points, but it happens in practice.  We also set S2Debug::DISABLE so
        // debug binaries won't abort on such polylines.
        testNearlyCovers("0:1, 0:2, 0:2, 0:3", "0:1, 0:1, 0:1, 0:3", 1e-10, true, true)
    }

    fun testS2PolylineCoveringTestCanChooseBetweenTwoPotentialStartingPoints() {
        // Can handle two possible starting points, only one of which leads to finding
        // a correct path.  In the first polyline, the edge from 0:1.1 to 0:0 and the
        // edge from 0:0.9 to 0:2 might be lucrative starting states for covering the
        // second polyline, because both edges are with the max_error of 1.5 degrees
        // from 0:10.  However, only the latter is actually effective.
        testNearlyCovers("0:11, 0:0, 0:9, 0:20", "0:10, 0:15", 1.5, false, true)
    }

    fun testS2PolylineCoveringTestStraightAndWigglyPolylinesCoverEachOther() {
        testNearlyCovers("40:1, 20:1",
                "39.9:0.9, 40:1.1, 30:1.15, 29:0.95, 28:1.1, 27:1.15, 26:1.05, 25:0.85, 24:1.1, 23:0.9, 20:0.99",
                0.2, true, true)
    }

    fun testS2PolylineCoveringTestMatchStartsAtLastVertex() {
        // The first polyline covers the second, but the matching segment starts at
        // the last vertex of the first polyline.
        testNearlyCovers(
                "0:0, 0:2", "0:2, 0:3", 1.5, false, true)
    }

    fun testS2PolylineCoveringTestMatchStartsAtDuplicatedLastVertex() {
        testNearlyCovers(
                "0:0, 0:2, 0:2, 0:2", "0:2, 0:3", 1.5, false, true)
    }

    fun testS2PolylineCoveringTestEmptyPolylines() {
        // We expect:
        //    anything.covers(empty) = true
        //    empty.covers(nonempty) = false
        testNearlyCovers("0:1, 0:2", "", 0.0, false, true)
        testNearlyCovers("", "", 0.0, true, true)
    }

}
