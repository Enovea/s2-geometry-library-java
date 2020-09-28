package dilivia.s2

import kotlin.math.nextTowards
import kotlin.math.pow

class S2EdgeCrosserTest : S2GeometryTestCase() {

// In non-debug builds, check that default-constructed and/or NaN S2Point
// arguments don't cause crashes, especially on the very first method call
// (since S2CopyingEdgeCrosser checks whether the first vertex of each edge is
// the same as the last vertex of the previous edged when deciding whether or
// not to call Restart).

    fun testCrossingSignInvalid(point: S2Point, expected: Int) {
        val crosser = S2EdgeCrosser(point, point)
        assertEquals(expected, crosser.crossingSign(point, point))
        val crosser2 = S2CopyingEdgeCrosser(point, point)
        assertEquals(expected, crosser2.crossingSign(point, point))
    }

    fun testEdgeOrVertexCrossingInvalid(point: S2Point, expected: Boolean) {
        val crosser = S2EdgeCrosser(point, point)
        assertEquals(expected, crosser.edgeOrVertexCrossing(point, point))
        val crosser2 = S2CopyingEdgeCrosser(point, point)
        assertEquals(expected, crosser2.edgeOrVertexCrossing(point, point))
    }

    /*
    fun testInvalidDefaultPoints() {
        // Check that default-constructed S2Point arguments don't cause crashes.
        val point = S2Point(0, 0, 0)
        testCrossingSignInvalid(point, 0)
        testEdgeOrVertexCrossingInvalid(point, false)
    }
    */

    /*
    fun testInvalidNanPoints() {
        // Check that NaN S2Point arguments don't cause crashes.
        val nan = Double.NaN
        val point = S2Point(nan, nan, nan)
        testCrossingSignInvalid(point, -1)
        testEdgeOrVertexCrossingInvalid(point, false)
    }
    */

    fun testCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point, robust: Int, edge_or_vertex: Boolean) {
        var robust = robust
        // Modify the expected result if two vertices from different edges match.
        if (a == c || a == d || b == c || b == d) robust = 0
        assertEquals(robust, S2EdgeCrossings.crossingSign(a, b, c, d))
        val crosser = S2EdgeCrosser(a, b, c)
        assertEquals(robust, crosser.crossingSign(d))
        assertEquals(robust, crosser.crossingSign(c))
        assertEquals(robust, crosser.crossingSign(d, c))
        assertEquals(robust, crosser.crossingSign(c, d))

        assertEquals(edge_or_vertex, S2EdgeCrossings.edgeOrVertexCrossing(a, b, c, d))
        crosser.restartAt(c)
        assertEquals(edge_or_vertex, crosser.edgeOrVertexCrossing(d))
        assertEquals(edge_or_vertex, crosser.edgeOrVertexCrossing(c))
        assertEquals(edge_or_vertex, crosser.edgeOrVertexCrossing(d, c))
        assertEquals(edge_or_vertex, crosser.edgeOrVertexCrossing(c, d))

        // Check that the crosser can be re-used.
        crosser.init(c, d)
        crosser.restartAt(a)
        assertEquals(robust, crosser.crossingSign(b))
        assertEquals(robust, crosser.crossingSign(a))

        // Now try all the same tests with CopyingEdgeCrosser.
        val crosser2 = S2CopyingEdgeCrosser(a, b, c)
        assertEquals(robust, crosser2.crossingSign(d))
        assertEquals(robust, crosser2.crossingSign(c))
        assertEquals(robust, crosser2.crossingSign(d, c))
        assertEquals(robust, crosser2.crossingSign(c, d))

        assertEquals(edge_or_vertex, S2EdgeCrossings.edgeOrVertexCrossing(a, b, c, d))
        crosser2.restartAt(c)
        assertEquals(edge_or_vertex, crosser2.edgeOrVertexCrossing(d))
        assertEquals(edge_or_vertex, crosser2.edgeOrVertexCrossing(c))
        assertEquals(edge_or_vertex, crosser2.edgeOrVertexCrossing(d, c))
        assertEquals(edge_or_vertex, crosser2.edgeOrVertexCrossing(c, d))

        // Check that the crosser can be re-used.
        crosser2.init(c, d)
        crosser2.restartAt(a)
        assertEquals(robust, crosser2.crossingSign(b))
        assertEquals(robust, crosser2.crossingSign(a))
    }

    fun testCrossings(a: S2Point, b: S2Point, c: S2Point, d: S2Point, robust: Int, edge_or_vertex: Boolean) {
        var a = a.normalize()
        var b = b.normalize()
        var c = c.normalize()
        var d = d.normalize()
        testCrossing(a, b, c, d, robust, edge_or_vertex)
        testCrossing(b, a, c, d, robust, edge_or_vertex)
        testCrossing(a, b, d, c, robust, edge_or_vertex)
        testCrossing(b, a, d, c, robust, edge_or_vertex)
        testCrossing(a, a, c, d, -1, false)
        testCrossing(a, b, c, c, -1, false)
        testCrossing(a, a, c, c, -1, false)
        testCrossing(a, b, a, b, 0, true)
        testCrossing(c, d, a, b, robust, edge_or_vertex != (robust == 0))
    }

    fun testCrossings() {
        // The real tests of edge crossings are in s2{loop,polygon}_test,
        // but we do a few simple tests here.

        // Two regular edges that cross.
        testCrossings(S2Point(1, 2, 1), S2Point(1.0, -3.0, 0.5), S2Point(1.0, -0.5, -3.0), S2Point(0.1, 0.5, 3.0), 1, true)

        // Two regular edges that intersect antipodal points.
        testCrossings(S2Point(1, 2, 1), S2Point(1.0, -3.0, 0.5), S2Point(-1.0, 0.5, 3.0), S2Point(-0.1, -0.5, -3.0), -1, false)

        // Two edges on the same great circle that start at antipodal points.
        testCrossings(S2Point(0, 0, -1), S2Point(0, 1, 0), S2Point(0, 0, 1), S2Point(0, 1, 1), -1, false)

        // Two edges that cross where one vertex is S2::Origin().
        testCrossings(S2Point(1, 0, 0), S2Point.origin(), S2Point(1.0, -0.1, 1.0), S2Point(1.0, 1.0, -0.1), 1, true)

        // Two edges that intersect antipodal points where one vertex is
        // S2::Origin().
        testCrossings(S2Point(1, 0, 0), S2Point.origin(), S2Point(-1.0, 0.1, -1.0), S2Point(-1.0, -1.0, 0.1), -1, false)

        // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
        // and edge CD is further CCW around (2,3,4) than AB.
        testCrossings(S2Point(2, 3, 4), S2Point(-1, 2, 5), S2Point(7, -2, 3), S2Point(2, 3, 4), 0, false)

        // Two edges that barely cross each other near the middle of one edge.  The
        // edge AB is approximately in the x=y plane, while CD is approximately
        // perpendicular to it and ends exactly at the x=y plane.
        testCrossings(S2Point(1, 1, 1), S2Point(1.0, 1.0.nextTowards(0.0), -1.0), S2Point(11, -12, -1), S2Point(10, 10, 1), 1, true)

        // In this version, the edges are separated by a distance of about 1e-15.
        testCrossings(S2Point(1, 1, 1), S2Point(1.0, 1.0.nextTowards(2.0), -1.0), S2Point(1, -1, 0), S2Point(1, 1, 0), -1, false)

        // Two edges that barely cross each other near the end of both edges.  This
        // example cannot be handled using regular double-precision arithmetic due
        // to floating-point underflow.
        testCrossings(S2Point(0, 0, 1), S2Point(2.0, -1e-323, 1.0), S2Point(1, -1, 1), S2Point(1e-323, 0.0, 1.0), 1, true)

        // In this version, the edges are separated by a distance of about 1e-640.
        testCrossings(S2Point(0, 0, 1), S2Point(2.0, 1e-323, 1.0), S2Point(1, -1, 1), S2Point(1e-323, 0.0, 1.0), -1, false)

        // Two edges that barely cross each other near the middle of one edge.
        // Computing the exact determinant of some of the triangles in this test
        // requires more than 2000 bits of precision.
        testCrossings(S2Point(1.0, -1e-323, -1e-323), S2Point(1e-323, 1.0, 1e-323), S2Point(1.0, -1.0, 1e-323), S2Point(1, 1, 0), 1, true)

        // In this version, the edges are separated by a distance of about 1e-640.
        testCrossings(S2Point(1.0, 1e-323, -1e-323), S2Point(-1e-323, 1.0, 1e-323), S2Point(1.0, -1.0, 1e-323), S2Point(1, 1, 0), -1, false)
    }

    fun testCollinearEdgesThatDontTouch() {
        val kIters = 500
        repeat(kIters) {
            val a = S2Random.randomPoint()
            val d = S2Random.randomPoint()
            val b = S2EdgeDistances.interpolate(0.05, a, d)
            val c = S2EdgeDistances.interpolate(0.95, a, d)
            assertTrue(0 >= S2EdgeCrossings.crossingSign(a, b, c, d))
            assertTrue(0 >= S2EdgeCrossings.crossingSign(a, b, c, d))
            val crosser = S2EdgeCrosser(a, b, c)
            assertTrue(0 >= crosser.crossingSign(d))
            assertTrue(0 >= crosser.crossingSign(c))
        }
    }

    fun testCoincidentZeroLengthEdgesThatDontTouch() {
        // It is important that the edge primitives can handle vertices that exactly
        // exactly proportional to each other, i.e. that are not identical but are
        // nevertheless exactly coincident when projected onto the unit sphere.
        // There are various ways that such points can arise.  For example,
        // Normalize() itself is not idempotent: there exist distinct points A,B
        // such that Normalize(A) == B  and Normalize(B) == A.  Another issue is
        // that sometimes calls to Normalize() are skipped when the result of a
        // calculation "should" be unit length mathematically (e.g., when computing
        // the cross product of two orthonormal vectors).
        //
        // This test checks pairs of edges AB and CD where A,B,C,D are exactly
        // coincident on the sphere and the norms of A,B,C,D are monotonically
        // increasing.  Such edge pairs should never intersect.  (This is not
        // obvious, since it depends on the particular symbolic perturbations used
        // by s2pred::Sign().  It would be better to replace this with a test that
        // says that the CCW results must be consistent with each other.)
        val kIters = 1000
        repeat(kIters) {
            // Construct a point P where every component is zero or a power of 2.
            val p = MutableS2Point()
            for (i in 0..2) {
                val binary_exp = S2Random.skewed(11)
                p[i] = if (binary_exp > 1022) 0.0 else 2.0.pow(-binary_exp.toDouble())
            }
            // If all components were zero, try again.  Note that normalization may
            // convert a non-zero point into a zero one due to underflow (!)
            if (p.norm2() == 0.0) return@repeat
            p.normalize()
            if (p[0] == 0.0 && p[1] == 0.0 && p[2] == 0.0) return@repeat

            // Now every non-zero component should have exactly the same mantissa.
            // This implies that if we scale the point by an arbitrary factor, every
            // non-zero component will still have the same mantissa.  Scale the points
            // so that they are all distinct and are still very likely to satisfy
            // S2::IsUnitLength (which allows for a small amount of error in the norm).
            val a = (1 - 3e-16) * p
            val b = (1 - 1e-16) * p
            val c = p
            val d = (1 + 2e-16) * p
            if (!S2Point.isUnitLength(a) || !S2Point.isUnitLength(d)) return@repeat
            // Verify that the expected edges do not cross.
            assertTrue(0 >= S2EdgeCrossings.crossingSign(a, b, c, d))
            val crosser = S2EdgeCrosser(a, b, c)
            assertTrue(0 >= crosser.crossingSign(d))
            assertTrue(0 >= crosser.crossingSign(c))
        }
    }
}

