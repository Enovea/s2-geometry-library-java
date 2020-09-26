package dilivia.s2

import com.google.common.geometry.S2.DBL_EPSILON
import dilivia.s2.math.R2Point
import mu.KotlinLogging
import java.lang.Math.pow
import kotlin.math.*

class S2EdgeClippingTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger {}

    fun testFaceClipping(a_raw: S2Point, b_raw: S2Point) {
        val a = a_raw.normalize()
        val b = b_raw.normalize()
        // TODO(ericv): Remove the following line once S2::RobustCrossProd is
        // extended to use simulation of simplicity.
        if (a == -b || S2Point.robustCrossProd(a, b).norm2() == 0.0) return

        if (!a.isUnitLength()) {
            logger.info { "testFaceClipping: a = $a is not unit length. |a|^2 = ${a.norm2()}" }
            return
        }

        if (!b.isUnitLength()) {
            logger.info { "testFaceClipping: b = $b is not unit length. |b|^2 = ${b.norm2()}" }
            return
        }

        // First we test GetFaceSegments.
        val segments = mutableListOf<S2EdgeClipping.FaceSegment>()
        S2EdgeClipping.getFaceSegments(a, b, segments)
        val n = segments.size
        assertTrue(n >= 1)

        if (logger.isTraceEnabled) {
            var msg = "\nA=$a_raw\nB=$b_raw"
            msg += "\nN=${S2Point.robustCrossProd(a, b)}\nSegments:\n"
            for ((i, s) in segments.withIndex()) {
                msg += "$i: face=${s.face}, a=${s.a}, b=${s.b}\n"
            }
            logger.trace { msg }
        }

        val biunit = R2Rect(R1Interval(-1, 1), R1Interval(-1, 1))
        val kErrorRadians = S2EdgeClipping.kFaceClipErrorRadians

        // The first and last vertices should approximately equal A and B.
        assertTrue(a.angle(S2Coords.faceUVtoXYZ(segments[0].face, segments[0].a)) <= kErrorRadians)
        assertTrue(b.angle(S2Coords.faceUVtoXYZ(segments[n - 1].face, segments[n - 1].b)) <= kErrorRadians)

        val norm = S2Point.robustCrossProd(a, b).normalize()
        val a_tangent = norm.crossProd(a)
        val b_tangent = b.crossProd(norm)
        segments.forEachIndexed { i, _ ->
            // Vertices may not protrude outside the biunit square.
            assertTrue(biunit.contains(segments[i].a))
            assertTrue(biunit.contains(segments[i].b))
            if (i != 0) {
                // The two representations of each interior vertex (on adjacent faces)
                // must correspond to exactly the same S2Point.
                assertTrue(segments[i - 1].face != segments[i].face)
                assertEquals(S2Coords.faceUVtoXYZ(segments[i - 1].face, segments[i - 1].b), S2Coords.faceUVtoXYZ(segments[i].face, segments[i].a))

                // Interior vertices should be in the plane containing A and B, and should
                // be contained in the wedge of angles between A and B (i.e., the dot
                // products with a_tangent and b_tangent should be non-negative).
                val p = S2Coords.faceUVtoXYZ(segments[i].face, segments[i].a).normalize()
                assertTrue(abs(p.dotProd(norm)) <= kErrorRadians)
                assertTrue(p.dotProd(a_tangent) >= -kErrorRadians)
                assertTrue(p.dotProd(b_tangent) >= -kErrorRadians)
            }
        }

        // Now we test ClipToPaddedFace (sometimes with a padding of zero).  We do
        // this by defining an (x,y) coordinate system for the plane containing AB,
        // and converting points along the great circle AB to angles in the range
        // [-Pi, Pi].  We then accumulate the angle intervals spanned by each
        // clipped edge; the union over all 6 faces should approximately equal the
        // interval covered by the original edge.
        val padding = if (S2Random.oneIn(10)) 0.0 else 1e-10 * pow(1e-5, S2Random.randomDouble())
        val x_axis = a
        val y_axis = a_tangent
        val expected_angles = S1Interval(0.0, a.angle(b))
        val max_angles = expected_angles.expanded(kErrorRadians)
        var actual_angles: S1Interval = S1Interval()
        for (face in 0..5) {
            val clippedAB = S2EdgeClipping.clipToPaddedFace(a, b, face, padding)
            if (clippedAB != null) {
              val (a_uv, b_uv) = clippedAB
                val a_clip = S2Coords.faceUVtoXYZ(face, a_uv).normalize()
                val b_clip = S2Coords.faceUVtoXYZ(face, b_uv).normalize()
                assertTrue(abs(a_clip.dotProd(norm)) <= kErrorRadians)
                assertTrue(abs(b_clip.dotProd(norm)) <= kErrorRadians)
                if (a_clip.angle(a) > kErrorRadians) {
                    assertEquals(1 + padding, max(abs(a_uv[0]), abs(a_uv[1])))
                }
                if (b_clip.angle(b) > kErrorRadians) {
                    assertEquals(1 + padding, max(abs(b_uv[0]), abs(b_uv[1])))
                }
                val a_angle = atan2(a_clip.dotProd(y_axis), a_clip.dotProd(x_axis))
                val b_angle = atan2(b_clip.dotProd(y_axis), b_clip.dotProd(x_axis))
                // Rounding errors may cause b_angle to be slightly less than a_angle.
                // We handle this by constructing the interval with FromPointPair(),
                // which is okay since the interval length is much less than M_PI.
                val face_angles = S1Interval.fromPointPair(a_angle, b_angle)
                assertTrue(max_angles.contains(face_angles))
                actual_angles = actual_angles.union(face_angles)
            }
        }
        assertTrue(actual_angles.expanded(kErrorRadians).contains(expected_angles))
    }

    fun testFaceClippingEdgePair(a: S2Point, b: S2Point) {
        testFaceClipping(a, b)
        testFaceClipping(b, a)
    }

    // This function is designed to choose line segment endpoints that are difficult
    // to handle correctly.  Given two adjacent cube vertices P and Q, it returns
    // either an edge midpoint, face midpoint, or corner vertex along the edge PQ
    // and then perturbs it slightly.  It also sometimes returns a random point from
    // anywhere on the sphere.
    fun perturbedCornerOrMidpoint(p: S2Point, q: S2Point): S2Point {
        var a = (S2Random.randomInt(3) - 1) * p + (S2Random.randomInt(3) - 1) * q
        if (S2Random.oneIn(10)) {
            // This perturbation often has no effect except on coordinates that are
            // zero, in which case the perturbed value is so small that operations on
            // it often result in underflow.
            a += 1e-300.pow(S2Random.randomDouble()) * S2Random.randomPoint()
        } else if (S2Random.oneIn(2)) {
            // For coordinates near 1 (say > 0.5), this perturbation yields values
            // that are only a few representable values away from the initial value.
            a += 4 * DBL_EPSILON * S2Random.randomPoint()
        } else {
            // A perturbation whose magnitude is in the range [1e-25, 1e-10].
            a += 1e-10 * 1e-15.pow(S2Random.randomDouble()) * S2Random.randomPoint()
        }
        if (a.norm2() < Double.MIN_VALUE) {
            // If a.Norm2() is denormalized, Normalize() loses too much precision.
            return perturbedCornerOrMidpoint(p, q)
        }
        return a
    }

    fun testFaceClipping() {
        // Start with a few simple cases.
        // An edge that is entirely contained within one cube face:
        testFaceClippingEdgePair(S2Point(1.0, -0.5, -0.5), S2Point(1.0, 0.5, 0.5))
        // An edge that crosses one cube edge:
        testFaceClippingEdgePair(S2Point(1, 0, 0), S2Point(0, 1, 0))
        // An edge that crosses two opposite edges of face 0:
        testFaceClippingEdgePair(S2Point(0.75, 0.0, -1.0), S2Point(0.75, 0.0, 1.0))
        // An edge that crosses two adjacent edges of face 2:
        testFaceClippingEdgePair(S2Point(1.0, 0.0, 0.75), S2Point(0.0, 1.0, 0.75))
        // An edge that crosses three cube edges (four faces):
        testFaceClippingEdgePair(S2Point(1.0, 0.9, 0.95), S2Point(-1.0, 0.95, 0.9))

        // Comprehensively test edges that are difficult to handle, especially those
        // that nearly follow one of the 12 cube edges.
        val biunit = R2Rect(R1Interval(-1, 1), R1Interval(-1, 1))
        val kIters = 1000;  // Test passes with 1e6 iterations
        repeat(kIters) { iter ->
            logger.trace { "Iteration $iter" }
            // Choose two adjacent cube corners P and Q.
            val face = S2Random.randomInt(6)
            val i = S2Random.randomInt(4)
            val j = (i + 1) and 3
            val p = S2Coords.faceUVtoXYZ(face, biunit.getVertex(i))
            val q = S2Coords.faceUVtoXYZ(face, biunit.getVertex(j))

            // Now choose two points that are nearly on the edge PQ, preferring points
            // that are near cube corners, face midpoints, or edge midpoints.
            val a = perturbedCornerOrMidpoint(p, q)
            val b = perturbedCornerOrMidpoint(p, q)
            testFaceClipping(a, b)
        }
    }

    // Choose a random point in the rectangle defined by points A and B, sometimes
    // returning a point on the edge AB or the points A and B themselves.
    fun chooseRectPoint(a: R2Point, b: R2Point): R2Point {
        return if (S2Random.oneIn(5)) {
            if (S2Random.oneIn(2)) a else b
        } else if (S2Random.oneIn(3)) {
            a + S2Random.randomDouble() * (b - a)
        } else {
            // a[i] may be >, <, or == b[i], so we write it like this instead
            // of using UniformDouble.
            R2Point(a[0] + S2Random.randomDouble() * (b[0] - a[0]),
                    a[1] + S2Random.randomDouble() * (b[1] - a[1]))
        }
    }

    // Given a point X on the line AB (which is checked), return the fraction "t"
    // such that x = (1-t)*a + t*b.  Return 0 if A = B.
    fun getFraction(x: R2Point, a: R2Point, b: R2Point): Double {
        // A bound for the error in edge clipping plus the error in the calculation
        // below (which is similar to IntersectsRect).
        val kError = (S2EdgeClipping.kEdgeClipErrorUVDist + S2EdgeClipping.kIntersectsRectErrorUVDist)
        if (a == b) return 0.0
        val dir = (b - a).normalize()
        assertTrue(abs((x - a).dotProd(dir.ortho())) <= kError)
        return (x - a).dotProd(dir)
    }

    // Given a point P representing a possibly clipped endpoint A of an edge AB,
    // verify that "clip" contains P, and that if clipping occurred (i.e., P != A)
    // then P is on the boundary of "clip".
    fun checkPointOnBoundary(p: R2Point, a: R2Point, clip: R2Rect) {
        assertTrue(clip.contains(p))
        if (p != a) {
            assertFalse(clip.contains(R2Point(p[0].nextTowards(a[0]), p[1].nextTowards(a[1]))))
        }
    }

    // Given an edge AB and a rectangle "clip", verify that IntersectsRect(),
    // ClipEdge(), and ClipEdgeBound() produce consistent results.
    fun testClipEdge(a: R2Point, b: R2Point, clip: R2Rect) {
        // A bound for the error in edge clipping plus the error in the
        // IntersectsRect calculation below.
        val kError = (S2EdgeClipping.kEdgeClipErrorUVDist + S2EdgeClipping.kIntersectsRectErrorUVDist)
        val clippedEdge = S2EdgeClipping.clipEdge(a, b, clip)
        if (clippedEdge == null) {
            assertFalse(S2EdgeClipping.intersectsRect(a, b, clip.expanded(-kError)))
        } else {
          val (a_clipped, b_clipped) = clippedEdge
            assertTrue(S2EdgeClipping.intersectsRect(a, b, clip.expanded(kError)))
            // Check that the clipped points lie on the edge AB, and that the points
            // have the expected order along the segment AB.
            assertTrue(getFraction(a_clipped, a, b) <= getFraction(b_clipped, a, b))
            // Check that the clipped portion of AB is as large as possible.
            checkPointOnBoundary(a_clipped, a, clip)
            checkPointOnBoundary(b_clipped, b, clip)
        }
        // Choose a random initial bound to pass to ClipEdgeBound.
        val initial_clip = R2Rect.fromPointPair(chooseRectPoint(a, b), chooseRectPoint(a, b))
        val bound = S2EdgeClipping.getClippedEdgeBound(a, b, initial_clip).toMutable()
        if (bound.isEmpty) return  // Precondition of ClipEdgeBound not met
        val max_bound = bound.intersection(clip)
        if (!S2EdgeClipping.clipEdgeBound(a, b, clip, bound)) {
            assertFalse(S2EdgeClipping.intersectsRect(a, b, max_bound.expanded(-kError)))
        } else {
            assertTrue(S2EdgeClipping.intersectsRect(a, b, max_bound.expanded(kError)))
            // Check that the bound is as large as possible.
            val ai = if (a[0] > b[0]) 1 else 0
            val aj = if (a[1] > b[1]) 1 else 0
            checkPointOnBoundary(bound.getVertex(ai, aj), a, max_bound)
            checkPointOnBoundary(bound.getVertex(1 - ai, 1 - aj), b, max_bound)
        }
    }

    // Given an interval "clip", randomly choose either a value in the interval, a
// value outside the interval, or one of the two interval endpoints, ensuring
// that all cases have reasonable probability for any interval "clip".
    fun chooseEndpoint(clip: R1Interval): Double {
        if (S2Random.oneIn(5)) {
            return if (S2Random.oneIn(2)) clip.lo else clip.hi
        } else {
            return when (S2Random.randomInt(3)) {
              0 -> clip.lo - S2Random.randomDouble()
              1 -> clip.hi + S2Random.randomDouble()
                else -> clip.lo + S2Random.randomDouble() * clip.length
            }
        }
    }

    // Given a rectangle "clip", choose a point that may lie in the rectangle
// interior, along an extended edge, exactly at a vertex, or in one of the
// eight regions exterior to "clip" that are separated by its extended edges.
// Also sometimes return points that are exactly on one of the extended
// diagonals of "clip".  All cases are reasonably likely to occur for any
// given rectangle "clip".
    fun chooseEndpoint(clip: R2Rect): R2Point {
        if (S2Random.oneIn(10)) {
            // Return a point on one of the two extended diagonals.
            val diag = S2Random.randomInt(2)
            val t = S2Random.randomDouble(-1.0, 2.0)
            return (1 - t) * clip.getVertex(diag) + t * clip.getVertex(diag + 2)
        } else {
            return R2Point(chooseEndpoint(clip[0]), chooseEndpoint(clip[1]))
        }
    }

    // Given a rectangle "clip", test the S2EdgeUtil edge clipping methods using
// many edges that are randomly constructed to trigger special cases.
    fun testEdgeClipping(clip: R2Rect) {
        val kIters = 1000;  // Test passes with 1e6 iterations
        repeat(kIters) { iter ->
            logger.trace { "Iteration $iter" }
            testClipEdge(chooseEndpoint(clip), chooseEndpoint(clip), clip)
        }
    }

    fun testEdgeClipping() {
        // Test clipping against random rectangles.
        for (i in 0..4) {
            testEdgeClipping(R2Rect.fromPointPair(
                    R2Point(
                            S2Random.randomDouble(-1.0, 1.0),
                            S2Random.randomDouble(-1.0, 1.0)
                    ),
                    R2Point(
                            S2Random.randomDouble(-1.0, 1.0),
                            S2Random.randomDouble(-1.0, 1.0)
                    )
            ))
        }
        // Also clip against one-dimensional, singleton, and empty rectangles.
        testEdgeClipping(R2Rect(R1Interval(-0.7, -0.7), R1Interval(0.3, 0.35)))
        testEdgeClipping(R2Rect(R1Interval(0.2, 0.5), R1Interval(0.3, 0.3)))
        testEdgeClipping(R2Rect(R1Interval(-0.7, 0.3), R1Interval(0, 0)))
        testEdgeClipping(R2Rect.fromPoint(R2Point(0.3, 0.8)))
        testEdgeClipping(R2Rect.empty())
    }

}