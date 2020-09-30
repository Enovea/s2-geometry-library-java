package dilivia.s2.region

import com.google.common.geometry.S2.M_PI
import com.google.common.geometry.S2.M_PI_2
import dilivia.s2.R1Interval
import dilivia.s2.S1Angle
import dilivia.s2.S2CellId
import dilivia.s2.S2EdgeCrossings
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2LatLng
import dilivia.s2.S2LatLngRect
import dilivia.s2.S2LatLngRectBounder
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.region.S2LoopTest.Companion.RelationFlags.CONTAINED
import dilivia.s2.region.S2LoopTest.Companion.RelationFlags.CONTAINS
import dilivia.s2.region.S2LoopTest.Companion.RelationFlags.COVERS
import dilivia.s2.region.S2LoopTest.Companion.RelationFlags.DISJOINT
import mu.KotlinLogging
import java.util.*
import kotlin.math.abs
import kotlin.math.acos
import kotlin.math.asin
import kotlin.math.cos
import kotlin.math.min
import kotlin.math.sin
import kotlin.math.tan


class S2LoopTest : S2GeometryTestCase() {

    companion object {

        private val logger = KotlinLogging.logger(S2LoopTest::class.java.name)

        // The set of all loops declared below.
        val allLoops = mutableListOf<S2Loop>()

        // Some standard loops to use in the tests (see descriptions below).
        // The empty loop.
        val empty: S2Loop = addLoop(S2Loop(S2Loop.kEmpty))

        // The full loop.
        val full: S2Loop = addLoop(S2Loop(S2Loop.kFull))

        // The northern hemisphere, defined using two pairs of antipodal points.
        val north_hemi: S2Loop = addLoop("0:-180, 0:-90, 0:0, 0:90")

        // The northern hemisphere, defined using three points 120 degrees apart.
        val north_hemi3: S2Loop = addLoop("0:-180, 0:-60, 0:60")

        // The southern hemisphere, defined using two pairs of antipodal points.
        val south_hemi: S2Loop = addLoop("0:90, 0:0, 0:-90, 0:-180")

        // The western hemisphere, defined using two pairs of antipodal points.
        val west_hemi: S2Loop = addLoop("0:-180, -90:0, 0:0, 90:0")

        // The eastern hemisphere, defined using two pairs of antipodal points.
        val east_hemi: S2Loop = addLoop("90:0, 0:0, -90:0, 0:-180")

        // The "near" hemisphere, defined using two pairs of antipodal points.
        val near_hemi: S2Loop = addLoop("0:-90, -90:0, 0:90, 90:0")

        // The "far" hemisphere, defined using two pairs of antipodal points.
        val far_hemi: S2Loop = addLoop("90:0, 0:90, -90:0, 0:-90")

        // A spiral stripe that slightly over-wraps the equator.
        val candy_cane: S2Loop = addLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")

        // A small clockwise loop in the northern & eastern hemisperes.
        val small_ne_cw: S2Loop = addLoop("35:20, 45:20, 40:25")

        // Loop around the north pole at 80 degrees.
        val arctic_80: S2Loop = addLoop("80:-150, 80:-30, 80:90")

        // Loop around the south pole at 80 degrees.
        val antarctic_80: S2Loop = addLoop("-80:120, -80:0, -80:-120")

        // A completely degenerate triangle along the equator that Sign()
        // considers to be CCW.
        val line_triangle: S2Loop = addLoop("0:1, 0:2, 0:3")

        // A nearly-degenerate CCW chevron near the equator with very long sides
        // (about 80 degrees).  Its area is less than 1e-640, which is too small
        // to represent in double precision.
        val skinny_chevron: S2Loop = addLoop("0:0, -1e-320:80, 0:1e-320, 1e-320:80")

        // A diamond-shaped loop around the point 0:180.
        val loop_a: S2Loop = addLoop("0:178, -1:180, 0:-179, 1:-180")

        // Another diamond-shaped loop around the point 0:180.
        val loop_b: S2Loop = addLoop("0:179, -1:180, 0:-178, 1:-180")

        // The intersection of A and B.
        val a_intersect_b: S2Loop = addLoop("0:179, -1:180, 0:-179, 1:-180")

        // The union of A and B.
        val a_union_b: S2Loop = addLoop("0:178, -1:180, 0:-178, 1:-180")

        // A minus B (concave).
        val a_minus_b: S2Loop = addLoop("0:178, -1:180, 0:179, 1:-180")

        // B minus A (concave).
        val b_minus_a: S2Loop = addLoop("0:-179, -1:180, 0:-178, 1:-180")

        // A shape gotten from A by adding a triangle to one edge, and
        // subtracting a triangle from the opposite edge.
        val loop_c: S2Loop = addLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")

        // A shape gotten from A by adding a triangle to one edge, and
        // adding another triangle to the opposite edge.
        val loop_d: S2Loop = addLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")

        //   3------------2
        //   |            |               ^
        //   |  7-8  b-c  |               |
        //   |  | |  | |  |      Latitude |
        //   0--6-9--a-d--1               |
        //   |  | |       |               |
        //   |  f-e       |               +----------->
        //   |            |                 Longitude
        //   4------------5
        //
        // Important: It is not okay to skip over collinear vertices when
        // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
        // uses symbolic perturbations to ensure that no three vertices are
        // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
        // collinear).  In other words, it is unpredictable (modulo knowing the
        // details of the symbolic perturbations) whether 0123 contains 06123,
        // for example.
        //
        // Loop E:  0,6,9,a,d,1,2,3
        // Loop F:  0,4,5,1,d,a,9,6
        // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
        // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
        // Loop I:  7,6,f,e,9,8
        val loop_e: S2Loop = addLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30")
        val loop_f: S2Loop = addLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34")
        val loop_g: S2Loop = addLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
        val loop_h: S2Loop = addLoop("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")
        val loop_i: S2Loop = addLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36")

        // Like loop_a, but the vertices are at leaf cell centers.
        val snapped_loop_a: S2Loop = addLoop(S2Loop(listOf<S2Point>(
                S2CellId.fromPoint(S2TextParser.makePoint("0:178")).toPoint(),
                S2CellId.fromPoint(S2TextParser.makePoint("-1:180")).toPoint(),
                S2CellId.fromPoint(S2TextParser.makePoint("0:-179")).toPoint(),
                S2CellId.fromPoint(S2TextParser.makePoint("1:-180")).toPoint()
        )))

        val kRectError = S2LatLngRectBounder.maxErrorForTests()


        fun addLoop(str: String): S2Loop {
            return addLoop(S2TextParser.makeLoop(str))
        }

        fun addLoop(loop: S2Loop): S2Loop {
            allLoops.add(loop)
            return loop
        }

        fun rotate(loop: S2Loop): S2Loop {
            val vertices = mutableListOf<S2Point>()
            for (i in 1..loop.numVertices()) {
                vertices.add(loop.vertex(i))
            }
            return S2Loop(vertices)
        }

        // Check that the curvature is *identical* when the vertex order is
        // rotated, and that the sign is inverted when the vertices are reversed.
        fun checkCurvatureInvariants(loop: S2Loop) {
            val expected = loop.getCurvature()
            var loop_copy = loop.clone()
            repeat(loop.numVertices()) {
                loop_copy.invert()
                assertEquals(-expected, loop_copy.getCurvature())
                loop_copy.invert()
                loop_copy = rotate(loop_copy)
                assertEquals(expected, loop_copy.getCurvature())
            }
        }


        // Checks that if a loop is normalized, it doesn't contain a
        // point outside of it, and vice versa.
        fun checkNormalizeAndContains(loop: S2Loop) {
            val p: S2Point = S2TextParser.makePoint("40:40")

            val flip = loop.clone()
            flip.invert()
            assertTrue(loop.isNormalized() xor loop.contains(p))
            assertTrue(flip.isNormalized() xor flip.contains(p))

            assertTrue(loop.isNormalized() xor flip.isNormalized())

            flip.normalize()
            assertFalse(flip.contains(p))
        }

        // Given a pair of loops where A contains B, check various identities.
        fun yestOneNestedPair(a: S2Loop, b: S2Loop) {
            assertTrue(a.contains(b))
            assertEquals(a.boundaryEquals(b), b.contains(a))
            assertEquals(!b.isEmpty(), a.intersects(b))
            assertEquals(!b.isEmpty(), b.intersects(a))
        }

        // Given a pair of disjoint loops A and B, check various identities.
        fun testOneDisjointPair(a: S2Loop, b: S2Loop) {
            assertFalse(a.intersects(b))
            assertFalse(b.intersects(a))
            assertEquals(b.isEmpty(), a.contains(b))
            assertEquals(a.isEmpty(), b.contains(a))
        }

        // Given loops A and B whose union covers the sphere, check various identities.
        fun testOneCoveringPair(a: S2Loop, b: S2Loop) {
            assertEquals(a.isFull(), a.contains(b))
            assertEquals(b.isFull(), b.contains(a))
            val a1 = (a.clone())
            a1.invert()
            val complementary = a1.boundaryEquals(b)
            assertEquals(!complementary, a.intersects(b))
            assertEquals(!complementary, b.intersects(a))
        }

        // Given loops A and B such that both A and its complement intersect both B
        // and its complement, check various identities.
        fun testOneOverlappingPair(a: S2Loop, b: S2Loop) {
            assertFalse(a.contains(b))
            assertFalse(b.contains(a))
            assertTrue(a.intersects(b))
            assertTrue(b.intersects(a))
        }

        // Given a pair of loops where A contains B, test various identities
        // involving A, B, and their complements.
        fun testNestedPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            val b1 = (b.clone())
            a1.invert()
            b1.invert()
            yestOneNestedPair(a, b)
            yestOneNestedPair(b1, a1)
            testOneDisjointPair(a1, b)
            testOneCoveringPair(a, b1)
        }

        // Given a pair of disjoint loops A and B, test various identities
        // involving A, B, and their complements.
        fun testDisjointPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            a1.invert()
            testNestedPair(a1, b)
        }

        // Given loops A and B whose union covers the sphere, test various identities
        // involving A, B, and their complements.
        fun testCoveringPair(a: S2Loop, b: S2Loop) {
            val b1 = (b.clone())
            b1.invert()
            testNestedPair(a, b1)
        }

        // Given loops A and B such that both A and its complement intersect both B
        // and its complement, test various identities involving these four loops.
        fun testOverlappingPair(a: S2Loop, b: S2Loop) {
            val a1 = (a.clone())
            val b1 = (b.clone())
            a1.invert()
            b1.invert()
            testOneOverlappingPair(a, b)
            testOneOverlappingPair(a1, b1)
            testOneOverlappingPair(a1, b)
            testOneOverlappingPair(a, b1)
        }

        object RelationFlags {
            val CONTAINS = 0x01  // A contains B
            val CONTAINED = 0x02  // B contains A
            val DISJOINT = 0x04  // A and B are disjoint (intersection is empty)
            val COVERS = 0x08  // (A union B) covers the entire sphere
        }

        // Verify the relationship between two loops A and B.  "flags" is the set of
        // RelationFlags that apply.  "shared_edge" means that the loops share at
        // least one edge (possibly reversed).
        fun testRelationWithDesc(a: S2Loop, b: S2Loop, flags: Int, shared_edge: Boolean, test_description: String) {
            logger.info { test_description }
            if (flags and CONTAINS != 0) {
                testNestedPair(a, b)
            }
            if (flags and CONTAINED != 0) {
                testNestedPair(b, a)
            }
            if (flags and COVERS != 0) {
                testCoveringPair(a, b)
            }
            if (flags and DISJOINT != 0) {
                testDisjointPair(a, b)
            } else if ((flags and (CONTAINS or CONTAINED or COVERS)) == 0) {
                testOverlappingPair(a, b)
            }
            if (!shared_edge && (flags and (CONTAINS or CONTAINED or DISJOINT)) != 0) {
                assertEquals(a.contains(b), a.containsNested(b))
            }
            // A contains the boundary of B if either A contains B, or the two loops
            // contain each other's boundaries and there are no shared edges (since at
            // least one such edge must be reversed, and therefore is not considered to
            // be contained according to the rules of CompareBoundary).
            var comparison = 0
            if ((flags and CONTAINS) != 0 || ((flags and COVERS) != 0 && !shared_edge)) {
                comparison = 1
            }
            // Similarly, A excludes the boundary of B if either A and B are disjoint,
            // or B contains A and there are no shared edges (since A is considered to
            // contain such edges according to the rules of CompareBoundary).
            if ((flags and DISJOINT) != 0 || ((flags and CONTAINED) != 0 && !shared_edge)) {
                comparison = -1
            }
            // CompareBoundary requires that neither loop is empty.
            if (!a.isEmpty() && !b.isEmpty()) {
                assertEquals(comparison, a.compareBoundary(b))
            }
        }

        fun testRelation(a: S2Loop, b: S2Loop, flags: Int, shared_edge: Boolean) = testRelationWithDesc(a, b, flags, shared_edge, "args $a , $b")

        fun makeCellLoop(begin: S2CellId, end: S2CellId): S2Loop {
            // Construct a CCW polygon whose boundary is the union of the cell ids
            // in the range [begin, end).  We add the edges one by one, removing
            // any edges that are already present in the opposite direction.
            val edges = TreeMap<S2Point, MutableSet<S2Point>>()
            var id = begin
            while (id != end) {
                val cell = S2Cell(id)
                for (k in 0..3) {
                    val a = cell.getVertex(k)
                    val b = cell.getVertex(k + 1)
                    val edgesB = edges[b]
                    val edgesA = edges[a]
                    if (edgesB == null || !edgesB.remove(a)) {
                        if (edgesA == null) edges[a] = mutableSetOf(b)
                        else edgesA.add(b)
                    } else if (edgesB.isEmpty()) {
                        edges.remove(b)
                    }
                }
                id = id.next()
            }

            // The remaining edges form a single loop.  We simply follow it starting
            // at an arbitrary vertex and build up a list of vertices.
            val vertices = mutableListOf<S2Point>()
            var p = edges.firstKey()
            while (!edges.isEmpty()) {
                assertEquals(1, edges.getValue(p).size)
                val next = edges.getValue(p).first()
                vertices.add(p)
                edges.remove(p)
                p = next
            }

            return S2Loop(vertices)
        }

        fun testNear(a_str: String, b_str: String, max_error: S1Angle, expected: Boolean) {
            val a = S2TextParser.makeLoop(a_str)
            val b = S2TextParser.makeLoop(b_str)
            assertEquals(a.boundaryNear(b, max_error), expected)
            assertEquals(b.boundaryNear(a, max_error), expected)
        }

        fun checkIdentical(loop: S2Loop, loop2: S2Loop) {
            assertEquals(loop.depth, loop2.depth)
            assertEquals(loop.numVertices(), loop2.numVertices())
            for (i in 0 until loop.numVertices()) {
                assertEquals(loop.vertex(i), loop2.vertex(i))
            }
            assertEquals(loop.isEmpty(), loop2.isEmpty())
            assertEquals(loop.isFull(), loop2.isFull())
            assertEquals(loop.depth, loop2.depth)
            assertEquals(loop.isNormalized(), loop2.isNormalized())
            assertEquals(loop.contains(S2Point.origin()), loop2.contains(S2Point.origin()))
            assertEquals(loop.rectBound, loop2.rectBound)
        }
    }

    fun testGetRectBound() {
        assertTrue(empty.rectBound.isEmpty)
        assertTrue(full.rectBound.isFull)
        assertTrue(candy_cane.rectBound.lng.isFull)
        assertTrue(candy_cane.rectBound.latLo().degrees() <= -20)
        assertTrue(candy_cane.rectBound.latHi().degrees() >= 10)
        assertTrue(small_ne_cw.rectBound.isFull)
        assertTrue(arctic_80.rectBound.approxEquals(S2LatLngRect(S2LatLng.fromDegrees(80, -180), S2LatLng.fromDegrees(90, 180)), kRectError))
        assertTrue(antarctic_80.rectBound.approxEquals(S2LatLngRect(S2LatLng.fromDegrees(-90, -180), S2LatLng.fromDegrees(-80, 180)), kRectError))

        // Create a loop that contains the complement of the "arctic_80" loop.
        val arctic_80inv = arctic_80.clone()
        arctic_80inv.invert()
        // The highest latitude of each edge is attained at its midpoint.
        val mid = (arctic_80inv.vertex(0) + arctic_80inv.vertex(1)) * 0.5
        assertEquals(arctic_80inv.rectBound.latHi().radians, S2LatLng.fromPoint(mid).lat().radians, kRectError.lat().radians)

        assertTrue(south_hemi.rectBound.lng.isFull)
        assertTrue(south_hemi.rectBound.lat.approxEquals(R1Interval(-M_PI_2, 0.0), kRectError.lat().radians))
    }

    fun testAreaConsistentWithCurvature() {
        // Check that the area computed using GetArea() is consistent with the
        // curvature of the loop computed using GetTurnAngle().  According to
        // the Gauss-Bonnet theorem, the area of the loop should be equal to 2*Pi
        // minus its curvature.
        for (loop in allLoops) {
            val area = loop.getArea()
            val gauss_area = 2 * M_PI - loop.getCurvature()
            // The error bound is sufficient for current tests but not guaranteed.
            assertTrue("Failed loop: ${loop.toString()}\nArea = $area, Gauss Area = $gauss_area",
                    abs(area - gauss_area) <= 1e-14
            )
        }
    }

    fun testGetAreaConsistentWithSign() {
        // Test that GetArea() returns an area near 0 for degenerate loops that
        // contain almost no points, and an area near 4*Pi for degenerate loops that
        // contain almost all points.
        val kMaxVertices = 6
        repeat(50) { i ->
            val num_vertices = 3 + S2Random.randomInt(kMaxVertices - 3 + 1)
            // Repeatedly choose N vertices that are exactly on the equator until we
            // find some that form a valid loop.
            var loop: S2Loop
            do {
                val vertices = mutableListOf<S2Point>()
                repeat(num_vertices) {
                    // We limit longitude to the range [0, 90] to ensure that the loop is
                    // degenerate (as opposed to following the entire equator).
                    vertices.add(S2LatLng.fromRadians(0.0, S2Random.randomDouble() * M_PI_2).toPoint())
                }
                loop = S2Loop(vertices, check = false)
            } while (!loop.isValid())
            val ccw = loop.isNormalized()
            assertEquals("Failed loop $i: $loop", if (ccw) 0.0 else 4 * M_PI, loop.getArea(), 1e-15)
            assertEquals(!ccw, loop.contains(S2Point(0, 0, 1)))
        }
    }

    fun testGetAreaAccuracy() {
        // TODO(ericv): Test that GetArea() has an accuracy significantly better
        // than 1e-15 on loops whose area is small.
    }

    fun testGetAreaAndCentroid() {
        assertEquals(0.0, empty.getArea())
        assertEquals(4 * M_PI, full.getArea())
        assertEquals(S2Point(0, 0, 0), empty.getCentroid())
        assertEquals(S2Point(0, 0, 0), full.getCentroid())

        assertEquals(north_hemi.getArea(), 2 * M_PI)
        assertEquals(east_hemi.getArea(), 2 * M_PI, 1e-15)

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices.  Then check that the area and centroid are
        // correct.

        repeat(50) {
            // Choose a coordinate frame for the spherical cap.
            val (x, y, z) = S2Random.randomFrame()

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want abs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            val kMaxDist = 1e-6
            val height = 2 * S2Random.randomDouble()
            val phi = asin(1 - height)
            var max_dtheta = 2 * acos(tan(abs(phi)) / tan(abs(phi) + kMaxDist))
            max_dtheta = min(M_PI, max_dtheta);  // At least 3 vertices.

            val vertices = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                vertices.add(x * (cos(theta) * cos(phi)) + y * (sin(theta) * cos(phi)) + z * sin(phi))
                theta += S2Random.randomDouble() * max_dtheta
            }
            val loop = S2Loop(vertices)
            val area = loop.getArea()
            val centroid = loop.getCentroid()
            val expected_area = 2 * M_PI * height
            assertTrue(abs(area - expected_area) <= 2 * M_PI * kMaxDist)
            val expected_centroid = z * (expected_area * (1 - 0.5 * height))
            assertTrue((centroid - expected_centroid).norm() <= 2 * kMaxDist)
        }
    }

    fun testGetCurvature() {
        assertEquals(2 * M_PI, empty.getCurvature())
        assertEquals(-2 * M_PI, full.getCurvature())

        assertEquals(0.0, north_hemi3.getCurvature(), 1e-15)
        checkCurvatureInvariants(north_hemi3)

        assertEquals(0.0, west_hemi.getCurvature(), 1e-15)
        checkCurvatureInvariants(west_hemi)

        // We don't have an easy way to estimate the curvature of this loop, but
        // we can still check that the expected invariants hold.
        checkCurvatureInvariants(candy_cane)

        assertEquals(2 * M_PI, line_triangle.getCurvature())
        checkCurvatureInvariants(line_triangle)

        assertEquals(2 * M_PI, skinny_chevron.getCurvature())
        checkCurvatureInvariants(skinny_chevron)

        // Build a narrow spiral loop starting at the north pole.  This is designed
        // to test that the error in GetCurvature is linear in the number of
        // vertices even when the partial sum of the curvatures gets very large.
        // The spiral consists of two "arms" defining opposite sides of the loop.
        val kArmPoints = 10000;    // Number of vertices in each "arm"
        val kArmRadius = 0.01;  // Radius of spiral.
        val vertices = Array<S2Point>(2 * kArmPoints) { S2Point() }
        vertices[kArmPoints] = S2Point(0, 0, 1)
        for (i in 0 until kArmPoints) {
            val angle = (2 * M_PI / 3) * i
            val x = cos(angle)
            val y = sin(angle)
            val r1 = i * kArmRadius / kArmPoints
            val r2 = (i + 1.5) * kArmRadius / kArmPoints
            vertices[kArmPoints - i - 1] = S2Point(r1 * x, r1 * y, 1.0).normalize()
            vertices[kArmPoints + i] = S2Point(r2 * x, r2 * y, 1.0).normalize()
        }
        // This is a pathological loop that contains many long parallel edges, and
        // takes tens of seconds to validate in debug mode.
        val spiral = S2Loop(vertices.toMutableList(), check = false)

        // Check that GetCurvature() is consistent with GetArea() to within the
        // error bound of the former.  We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that.  The error in GetArea() can be ignored for the
        // purposes of this test since it is generally much smaller.
        assertEquals(2 * M_PI - spiral.getArea(), spiral.getCurvature(), 0.01 * spiral.getCurvatureMaxError())
    }

    fun testNormalizedCompatibleWithContains() {
        checkNormalizeAndContains(line_triangle)
        checkNormalizeAndContains(skinny_chevron)
    }

    fun testContains() {
        // Check the full and empty loops have the correct containment relationship
        // with the special "vertex" that defines them.
        assertFalse(empty.contains(S2Loop.kEmpty[0]))
        assertTrue(full.contains(S2Loop.kFull[0]))

        assertTrue(candy_cane.contains(S2LatLng.fromDegrees(5, 71).toPoint()))

        // Create copies of these loops so that we can change the vertex order.
        var north_copy = (north_hemi.clone())
        var south_copy = (south_hemi.clone())
        var west_copy = (west_hemi.clone())
        var east_copy = (east_hemi.clone())
        repeat(4) {
            assertTrue(north_copy.contains(S2Point(0, 0, 1)))
            assertFalse(north_copy.contains(S2Point(0, 0, -1)))
            assertFalse(south_copy.contains(S2Point(0, 0, 1)))
            assertTrue(south_copy.contains(S2Point(0, 0, -1)))
            assertFalse(west_copy.contains(S2Point(0, 1, 0)))
            assertTrue(west_copy.contains(S2Point(0, -1, 0)))
            assertTrue(east_copy.contains(S2Point(0, 1, 0)))
            assertFalse(east_copy.contains(S2Point(0, -1, 0)))
            north_copy = rotate(north_copy)
            south_copy = rotate(south_copy)
            east_copy = rotate(east_copy)
            west_copy = rotate(west_copy)
        }

        // This code checks each cell vertex is contained by exactly one of
        // the adjacent cells.
        for (level in 0..2) {
            val loops = mutableListOf<S2Loop>()
            val loop_vertices = mutableListOf<S2Point>()
            val points = mutableSetOf<S2Point>()
            var id = S2CellId.begin(level)
            while (id != S2CellId.end(level)) {
                val cell = S2Cell(id)
                points.add(cell.getCenter())
                for (k in 0..3) {
                    loop_vertices.add(cell.getVertex(k))
                    points.add(cell.getVertex(k))
                }
                loops.add(S2Loop(loop_vertices))
                loop_vertices.clear()
                id = id.next()
            }
            for (point in points) {
                var count = 0
                for (loop in loops) {
                    if (loop.contains(point)) ++count
                }
                assertEquals(count, 1)
            }
        }
    }

    fun testContainsMatchesCrossingSign() {
        // This test demonstrates a former incompatibility between CrossingSign()
        // and Contains(const S2Point&).  It constructs an S2Cell-based loop L and
        // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
        // previously, Contains() returned false for both endpoints of E.
        //
        // The reason for the bug was that the loop bound was sometimes too tight.
        // The Contains() code for a0 bailed out early because a0 was found not to
        // be inside the bound of L.

        // Start with a cell that ends up producing the problem.
        val cell_id = S2CellId.fromPoint(S2Point(1, 1, 1)).parent(21)

        val children = Array(4) { S2Cell() }
        S2Cell(cell_id).subdivide(children)

        val points = Array<S2Point>(4) { S2Point() }
        for (i in 0..3) {
            // Note extra normalization. GetCenter() is already normalized.
            // The test results will no longer be inconsistent if the extra
            // Normalize() is removed.
            points[i] = children[i].getCenter().normalize()
        }

        val loop = S2Loop(points.toMutableList())

        // Get a vertex from a grandchild cell.
        // +---------------+---------------+
        // |               |               |
        // |    points[3]  |   points[2]   |
        // |       v       |       v       |
        // |       +-------+------ +       |
        // |       |       |       |       |
        // |       |       |       |       |
        // |       |       |       |       |
        // +-------+-------+-------+-------+
        // |       |       |       |       |
        // |       |    <----------------------- grandchild_cell
        // |       |       |       |       |
        // |       +-------+------ +       |
        // |       ^       |       ^       | <-- cell
        // | points[0]/a0  |     points[1] |
        // |               |               |
        // +---------------+---------------+
        val grandchild_cell = S2Cell(cell_id.child(0).child(2))
        val a0 = grandchild_cell.getVertex(0)

        // If this doesn't hold, the rest of the test is pointless.
        assertTrue(
                "This test depends on rounding errors that should make a0 slightly different from points[0]\npoints[0]: ${points[0]}\n       a0: $a0",
                points[0] != a0)

        // The edge from a0 to the origin crosses one boundary.
        assertEquals(-1, S2EdgeCrossings.crossingSign(a0, S2Point.origin(), loop.vertex(0), loop.vertex(1)))
        assertEquals(1, S2EdgeCrossings.crossingSign(a0, S2Point.origin(), loop.vertex(1), loop.vertex(2)))
        assertEquals(-1, S2EdgeCrossings.crossingSign(a0, S2Point.origin(), loop.vertex(2), loop.vertex(3)))
        assertEquals(-1, S2EdgeCrossings.crossingSign(a0, S2Point.origin(), loop.vertex(3), loop.vertex(4)))

        // Contains should return false for the origin, and true for a0.
        assertFalse(loop.contains(S2Point.origin()))
        assertTrue(loop.contains(a0))

        // Since a0 is inside the loop, it should be inside the bound.
        val bound = loop.rectBound
        assertTrue(bound.contains(a0))
    }


    fun testLoopRelations() {
        // Check full and empty relationships with normal loops and each other.
        testRelation(full, full, CONTAINS or CONTAINED or COVERS, true)
        testRelation(full, north_hemi, CONTAINS or COVERS, false)
        testRelation(full, empty, CONTAINS or DISJOINT or COVERS, false)
        testRelation(north_hemi, full, CONTAINED or COVERS, false)
        testRelation(north_hemi, empty, CONTAINS or DISJOINT, false)
        testRelation(empty, full, CONTAINED or DISJOINT or COVERS, false)
        testRelation(empty, north_hemi, CONTAINED or DISJOINT, false)
        testRelation(empty, empty, CONTAINS or CONTAINED or DISJOINT, false)

        testRelation(north_hemi, north_hemi, CONTAINS or CONTAINED, true)
        testRelation(north_hemi, south_hemi, DISJOINT or COVERS, true)
        testRelation(north_hemi, east_hemi, 0, false)
        testRelation(north_hemi, arctic_80, CONTAINS, false)
        testRelation(north_hemi, antarctic_80, DISJOINT, false)
        testRelation(north_hemi, candy_cane, 0, false)

        // We can't compare north_hemi3 vs. north_hemi or south_hemi because the
        // result depends on the "simulation of simplicity" implementation details.
        testRelation(north_hemi3, north_hemi3, CONTAINS or CONTAINED, true)
        testRelation(north_hemi3, east_hemi, 0, false)
        testRelation(north_hemi3, arctic_80, CONTAINS, false)
        testRelation(north_hemi3, antarctic_80, DISJOINT, false)
        testRelation(north_hemi3, candy_cane, 0, false)

        testRelation(south_hemi, north_hemi, DISJOINT or COVERS, true)
        testRelation(south_hemi, south_hemi, CONTAINS or CONTAINED, true)
        testRelation(south_hemi, far_hemi, 0, false)
        testRelation(south_hemi, arctic_80, DISJOINT, false)
        testRelation(south_hemi, antarctic_80, CONTAINS, false)
        testRelation(south_hemi, candy_cane, 0, false)

        testRelation(candy_cane, north_hemi, 0, false)
        testRelation(candy_cane, south_hemi, 0, false)
        testRelation(candy_cane, arctic_80, DISJOINT, false)
        testRelation(candy_cane, antarctic_80, DISJOINT, false)
        testRelation(candy_cane, candy_cane, CONTAINS or CONTAINED, true)

        testRelation(near_hemi, west_hemi, 0, false)

        testRelation(small_ne_cw, south_hemi, CONTAINS, false)
        testRelation(small_ne_cw, west_hemi, CONTAINS, false)

        testRelation(small_ne_cw, north_hemi, COVERS, false)
        testRelation(small_ne_cw, east_hemi, COVERS, false)

        testRelation(loop_a, loop_a, CONTAINS or CONTAINED, true)
        testRelation(loop_a, loop_b, 0, false)
        testRelation(loop_a, a_intersect_b, CONTAINS, true)
        testRelation(loop_a, a_union_b, CONTAINED, true)
        testRelation(loop_a, a_minus_b, CONTAINS, true)
        testRelation(loop_a, b_minus_a, DISJOINT, true)

        testRelation(loop_b, loop_a, 0, false)
        testRelation(loop_b, loop_b, CONTAINS or CONTAINED, true)
        testRelation(loop_b, a_intersect_b, CONTAINS, true)
        testRelation(loop_b, a_union_b, CONTAINED, true)
        testRelation(loop_b, a_minus_b, DISJOINT, true)
        testRelation(loop_b, b_minus_a, CONTAINS, true)

        testRelation(a_intersect_b, loop_a, CONTAINED, true)
        testRelation(a_intersect_b, loop_b, CONTAINED, true)
        testRelation(a_intersect_b, a_intersect_b, CONTAINS or CONTAINED, true)
        testRelation(a_intersect_b, a_union_b, CONTAINED, false)
        testRelation(a_intersect_b, a_minus_b, DISJOINT, true)
        testRelation(a_intersect_b, b_minus_a, DISJOINT, true)

        testRelation(a_union_b, loop_a, CONTAINS, true)
        testRelation(a_union_b, loop_b, CONTAINS, true)
        testRelation(a_union_b, a_intersect_b, CONTAINS, false)
        testRelation(a_union_b, a_union_b, CONTAINS or CONTAINED, true)
        testRelation(a_union_b, a_minus_b, CONTAINS, true)
        testRelation(a_union_b, b_minus_a, CONTAINS, true)

        testRelation(a_minus_b, loop_a, CONTAINED, true)
        testRelation(a_minus_b, loop_b, DISJOINT, true)
        testRelation(a_minus_b, a_intersect_b, DISJOINT, true)
        testRelation(a_minus_b, a_union_b, CONTAINED, true)
        testRelation(a_minus_b, a_minus_b, CONTAINS or CONTAINED, true)
        testRelation(a_minus_b, b_minus_a, DISJOINT, false)

        testRelation(b_minus_a, loop_a, DISJOINT, true)
        testRelation(b_minus_a, loop_b, CONTAINED, true)
        testRelation(b_minus_a, a_intersect_b, DISJOINT, true)
        testRelation(b_minus_a, a_union_b, CONTAINED, true)
        testRelation(b_minus_a, a_minus_b, DISJOINT, false)
        testRelation(b_minus_a, b_minus_a, CONTAINS or CONTAINED, true)
    }

    // Make sure the relations are correct if the loop crossing happens on
    // two ends of a shared boundary segment.
    fun testLoopRelationsWhenSameExceptPiecesStickingOutAndIn() {
        testRelation(loop_a, loop_c, 0, true)
        testRelation(loop_c, loop_a, 0, true)
        testRelation(loop_a, loop_d, CONTAINED, true)
        testRelation(loop_d, loop_a, CONTAINS, true)
        testRelation(loop_e, loop_f, DISJOINT, true)
        testRelation(loop_e, loop_g, CONTAINS, true)
        testRelation(loop_e, loop_h, 0, true)
        testRelation(loop_e, loop_i, 0, false)
        testRelation(loop_f, loop_g, DISJOINT, true)
        testRelation(loop_f, loop_h, 0, true)
        testRelation(loop_f, loop_i, 0, false)
        testRelation(loop_g, loop_h, CONTAINED, true)
        testRelation(loop_h, loop_g, CONTAINS, true)
        testRelation(loop_g, loop_i, DISJOINT, true)
        testRelation(loop_h, loop_i, CONTAINS, true)
    }


    fun testLoopRelations2() {
        // Construct polygons consisting of a sequence of adjacent cell ids
        // at some fixed level.  Comparing two polygons at the same level
        // ensures that there are no T-vertices.
        repeat(1000) {
            var begin = S2CellId((S2Random.randomLong() or 1L).toULong())
            if (!begin.isValid()) return@repeat
            begin = begin.parent(S2Random.randomInt(S2CellId.kMaxLevel))
            val a_begin = begin.advance(S2Random.skewed(6).toLong())
            val a_end = a_begin.advance(S2Random.skewed(6).toLong() + 1L)
            val b_begin = begin.advance(S2Random.skewed(6).toLong())
            val b_end = b_begin.advance(S2Random.skewed(6).toLong() + 1L)
            if (!a_end.isValid() || !b_end.isValid()) return@repeat

            val a = (makeCellLoop(a_begin, a_end))
            val b = (makeCellLoop(b_begin, b_end))
            val contained = (a_begin <= b_begin && b_end <= a_end)
            val intersects = (a_begin < b_end && b_begin < a_end)
            logger.trace { "Checking ${a.numVertices()} vs. ${b.numVertices()}, contained = $contained, intersects = $intersects" }
            assertEquals(a.contains(b), contained)
            assertEquals(a.intersects(b), intersects)
        }
    }

    fun testBoundsForLoopContainment() {
        // To reliably test whether one loop contains another, the bounds of the
        // outer loop are expanded slightly.  This test constructs examples where
        // this expansion is necessary and verifies that it is sufficient.
        repeat(1000) {
            // We construct a triangle ABC such that A,B,C are nearly colinear, B is
            // the point of maximum latitude, and the edge AC passes very slightly
            // below B (i.e., ABC is CCW).
            val b = (S2Random.randomPoint() + S2Point(0, 0, 1)).normalize()
            val v = b.crossProd(S2Point(0, 0, 1)).normalize()
            val a = S2EdgeDistances.interpolate(S2Random.randomDouble(), -v, b)
            val c = S2EdgeDistances.interpolate(S2Random.randomDouble(), b, v)
            if (S2Predicates.sign(a, b, c) < 0) return@repeat
            // Now construct another point D directly below B, and create two loops
            // ABCD and ACD.
            val d = S2Point(b.x(), b.y(), 0.0).normalize()
            val vertices = arrayOf(c, d, a, b)  // Reordered for convenience
            val outer = S2Loop(vertices.toMutableList())
            val inner = S2Loop(vertices.slice(0..2).toMutableList())
            // Now because the bounds calculation is less accurate when the maximum is
            // attained along an edge (rather than at a vertex), sometimes the inner
            // loop will have a *larger* bounding box than the outer loop.  We look
            // only for those cases.
            if (outer.rectBound.contains(inner.rectBound)) return@repeat
            assertTrue(outer.contains(inner))
        }
    }

    fun debugDumpCrossings(loop: S2Loop) {
        // This function is useful for debugging.
        println("""
            |Ortho(v1): ${S2Point.ortho(loop.vertex(1))}
            |Contains(kOrigin): ${loop.contains(S2Point.origin())}
        """.trimMargin())
        for (i in 1..loop.numVertices()) {
            val a = S2Point.ortho(loop.vertex(i))
            val b = loop.vertex(i - 1)
            val c = loop.vertex(i + 1)
            val o = loop.vertex(i)

            println("Vertex %d: [%.17g, %.17g, %.17g], %d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %d".format(
                    i, loop.vertex(i).x(), loop.vertex(i).y(), loop.vertex(i).z(),
                    i - 1, i, S2Predicates.sign(b, o, a),
                    i + 1, i, i - 1, S2Predicates.sign(c, o, b),
                    i, i + 1, S2Predicates.sign(a, o, c),
                    S2Predicates.orderedCCW(a, b, c, o)
            ))
        }
        for (i in 0 until loop.numVertices() + 2) {
            var orig = S2Point.origin()
            val dest: S2Point
            if (i < loop.numVertices()) {
                dest = loop.vertex(i)
                println("Origin.%d crosses:".format(i))
            } else {
                dest = S2Point(0, 0, 1)
                if (i == loop.numVertices() + 1) orig = loop.vertex(1)
                println("Case %d:".format(i))
            }
            for (j in 0 until loop.numVertices()) {
                println(" %d".format(S2EdgeCrossings.edgeOrVertexCrossing(orig, dest, loop.vertex(j), loop.vertex(j + 1))))
            }
            println()
        }
        for (i in 0..2 step 2) {
            println("Origin.v1 crossing v%d.v1: ".format(i))
            val a = S2Point.ortho(loop.vertex(1))
            val b = loop.vertex(i)
            val c = S2Point.origin()
            val o = loop.vertex(1)
            println("%d1R=%d, M1%d=%d, R1M=%d, crosses: %d".format(
                    i, S2Predicates.sign(b, o, a), i, S2Predicates.sign(c, o, b), S2Predicates.sign(a, o, c),
                    S2EdgeCrossings.edgeOrVertexCrossing(c, o, b, a))
            )
        }
    }

    fun testBoundaryNear() {
        val degree = S1Angle.degrees(1)

        testNear("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", degree * 0.5, true)
        testNear("0:0, 0:3, 0:7, 0:10, 3:7, 5:5", "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1", S1Angle.radians(1e-3), true)

        // All vertices close to some edge, but not equivalent.
        testNear("0:0, 0:2, 2:2, 2:0", "0:0, 1.9999:1, 0:2, 2:2, 2:0", degree * 0.5, false)

        // Two triangles that backtrack a bit on different edges.  A simple
        // greedy matching algorithm would fail on this example.
        val t1 = "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, 2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2"
        val t2 = "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, 0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4"
        testNear(t1, t2, degree * 1.5, true)
        testNear(t1, t2, degree * 0.5, false)
    }
/*

    fun testEncodeDecode(const S2Loop& loop)
    {
        Encoder encoder
                loop.Encode(& encoder)
        Decoder decoder (encoder.base(), encoder.length())
        S2Loop loop2
                loop2.set_s2debug_override(loop.s2debug_override())
        ASSERT_TRUE(loop2.Decode(& decoder))
        CheckIdentical(loop, loop2)
    }

    fun testEncodeDecode()
    {
        unique_ptr<S2Loop> l (s2textformat::MakeLoop("30:20, 40:20, 39:43, 33:35"))
        l.set_depth(3)
        TestEncodeDecode(*l)

        S2Loop empty (S2Loop.kEmpty)
        TestEncodeDecode(empty)
        S2Loop full (S2Loop.kFull)
        TestEncodeDecode(full)

        S2Loop uninitialized
                TestEncodeDecode(uninitialized)
    }
    */

    /*

    fun testEmptyFullSnapped(const S2Loop& loop, int level)
    {
        S2_CHECK(loop.is_emptyor_full())
        S2CellId cellid = S2CellId (loop.vertex(0)).parent(level)
        vector<S2Point> vertices = { cellid.toPoint() }
        S2Loop loop2 (vertices)
        assertTrue(loop.BoundaryEquals(& loop2))
        assertTrue(loop.BoundaryapproxEquals(loop2))
        assertTrue(loop.BoundaryNear(loop2))
    }

// Test converting the empty/full loops to S2LatLng representations.  (We
// don't bother testing E5/E6/E7 because that test is less demanding.)
    fun testEmptyFullLatLng(const S2Loop& loop)
    {
        S2_CHECK(loop.is_emptyor_full())
        vector<S2Point> vertices = { S2LatLng(loop.vertex(0)).toPoint() }
        S2Loop loop2 (vertices)
        assertTrue(loop.BoundaryEquals(& loop2))
        assertTrue(loop.BoundaryapproxEquals(loop2))
        assertTrue(loop.BoundaryNear(loop2))
    }

    fun testEmptyFullConversions(const S2Loop& loop)
    {
        TestEmptyFullSnapped(loop, S2CellId::kMaxLevel)
        TestEmptyFullSnapped(loop, 1);  // Worst case for approximation
        TestEmptyFullSnapped(loop, 0)
        TestEmptyFullLatLng(loop)
    }

    fun testEmptyFullLossyConversions()
    {
        // Verify that the empty and full loops can be encoded lossily.
        S2Loop empty (S2Loop.kEmpty)
        TestEmptyFullConversions(empty)

        S2Loop full (S2Loop.kFull)
        TestEmptyFullConversions(full)
    }

    fun testEncodeDecodeWithinScope()
    {
        unique_ptr<S2Loop> l (s2textformat::MakeLoop("30:20, 40:20, 39:43, 33:35"))
        l.set_depth(3)
        Encoder encoder
                l.Encode(& encoder)
        Decoder decoder1 (encoder.base(), encoder.length())

        // Initialize the loop using DecodeWithinScope and check that it is the
        // same as the original loop.
        S2Loop loop1
                ASSERT_TRUE(loop1.DecodeWithinScope(& decoder1))
        assertTrue(l.BoundaryEquals(& loop1))
        assertEquals(l.depth(), loop1.depth())
        assertEquals(l.rectBound, loop1.rectBound)

        // Initialize the same loop using Init with a vector of vertices, and
        // check that it doesn't deallocate the original memory.
        vector<S2Point> vertices = {
            loop1.vertex(0), loop1.vertex(2),
            loop1.vertex(3)
        }
        loop1.Init(vertices)
        Decoder decoder2 (encoder.base(), encoder.length())
        S2Loop loop2
                ASSERT_TRUE(loop2.DecodeWithinScope(& decoder2))
        assertTrue(l.BoundaryEquals(& loop2))
        assertEquals(l.vertex(1), loop2.vertex(1))
        EXPECT_NE(loop1.vertex(1), loop2.vertex(1))

        // Initialize loop2 using Decode with a decoder on different data.
        // Check that the original memory is not deallocated or overwritten.
        unique_ptr<S2Loop> l2 (s2textformat::MakeLoop("30:40, 40:75, 39:43, 80:35"))
        l2.set_depth(2)
        Encoder encoder2
                l2.Encode(& encoder2)
        Decoder decoder3 (encoder2.base(), encoder2.length())
        ASSERT_TRUE(loop2.Decode(& decoder3))
        Decoder decoder4 (encoder.base(), encoder.length())
        ASSERT_TRUE(loop1.DecodeWithinScope(& decoder4))
        assertTrue(l.BoundaryEquals(& loop1))
        assertEquals(l.vertex(1), loop1.vertex(1))
        EXPECT_NE(loop1.vertex(1), loop2.vertex(1))
    }

    fun testFourVertexCompressedLoopRequires36Bytes() {
        Encoder encoder
                TestEncodeCompressed(*snapped_loop_a, S2CellId::kMaxLevel, & encoder)

        // 1 byte for num_vertices
        // 1 byte for origin_inside and boolean indicating we did not
        //   encode the bound
        // 1 byte for depth
        // Vertices:
        // 1 byte for faces
        // 8 bytes for each vertex.
        // 1 byte indicating that there is no unsnapped vertex.
        assertEquals(37, encoder.length())
    }

    fun testCompressedEncodedLoopDecodesApproxEqual() {
        unique_ptr<S2Loop> loop (snapped_loop_a.clone())
        loop.set_depth(3)

        Encoder encoder
                TestEncodeCompressed(*loop, S2CellId::kMaxLevel, & encoder)
        S2Loop decoded_loop
                TestDecodeCompressed(encoder, S2CellId::kMaxLevel, & decoded_loop)
        CheckIdentical(*loop, decoded_loop)
    }

// This test checks that S2Loops created directly from S2Cells behave
// identically to S2Loops created from the vertices of those cells; this
// previously was not the case, because S2Cells calculate their bounding
// rectangles slightly differently, and S2Loops created from them just copied
// the S2Cell bounds.
    fun testS2CellConstructorAndContains()
    {
        S2Cell cell (S2CellId(S2LatLng::FromE6(40565459, -74645276)))
        S2Loop cell_as_loop (cell)

        vector<S2Point> vertices
                for (int i = 0; i < cell_as_loop.num_vertices(); ++i) {
        vertices.push_back(cell_as_loop.vertex(i))
    }
        S2Loop loop_copy (vertices)
        assertTrue(loop_copy.contains(& cell_as_loop))
        assertTrue(cell_as_loop.contains(& loop_copy))

        // Demonstrates the reason for this test; the cell bounds are more
        // conservative than the resulting loop bounds.
        assertFalse(loop_copy.rectBound.contains(cell.rectBound))
    }

// Construct a loop using s2textformat::MakeLoop(str) and check that it
// produces a validation error that includes "snippet".
    static void CheckLoopIsInvalid(const string& str, const string& snippet)
    {
        unique_ptr<S2Loop> loop (s2textformat::MakeLoop(str, S2Debug::DISABLE))
        S2Error error
                assertTrue(loop.FindValidationError(& error))
        EXPECT_NE(string::npos, error.text().find(snippet))
    }

    static void CheckLoopIsInvalid(const vector<S2Point>& points,
    const string& snippet)
    {
        S2Loop l (points, S2Debug::DISABLE)
        S2Error error
                assertTrue(l.FindValidationError(& error))
        EXPECT_NE(string::npos, error.text().find(snippet))
    }

    fun testIsValidDetectsInvalidLoops()
    {
        // Not enough vertices.  Note that all single-vertex loops are valid; they
        // are interpreted as being either empty or full.
        CheckLoopIsInvalid("", "at least 3 vertices")
        CheckLoopIsInvalid("20:20, 21:21", "at least 3 vertices")

        // There is a degenerate edge
        CheckLoopIsInvalid("20:20, 20:20, 20:21", "degenerate")
        CheckLoopIsInvalid("20:20, 20:21, 20:20", "degenerate")

        // There is a duplicate vertex
        CheckLoopIsInvalid("20:20, 21:21, 21:20, 20:20, 20:21", "duplicate vertex")

        // Some edges cross
        CheckLoopIsInvalid("20:20, 21:21, 21:20.5, 21:20, 20:21", "crosses")

        // Points with non-unit length (triggers S2_DCHECK failure in debug)
        EXPECT_DEBUG_DEATH(
                CheckLoopIsInvalid({ S2Point(2, 0, 0), S2Point(0, 1, 0), S2Point(0, 0, 1) },
                        "unit length"),
                "IsUnitLength")

        // Adjacent antipodal vertices
        CheckLoopIsInvalid({ S2Point(1, 0, 0), S2Point(-1, 0, 0), S2Point(0, 0, 1) },
                "antipodal")
    }

// Helper function for testing the distance methods.  "boundary_x" is the
// expected result of projecting "x" onto the loop boundary.  For convenience
// it can be set to S2Point() to indicate that (boundary_x == x).
    fun testDistanceMethods(const S2Loop& loop, const S2Point& x,
    S2Point boundary_x)
    {
        // This error is not guaranteed by the implementation but is okay for tests.
        const S1Angle kMaxError = S1Angle::Radians(1e-15)

        if (boundary_x == S2Point()) boundary_x = x
        EXPECT_LE(S1Angle(boundary_x, loop.ProjectToBoundary(x)), kMaxError)

        if (loop.is_emptyor_full()) {
            assertEquals(S1Angle::Infinity(), loop.GetDistanceToBoundary(x))
        } else {
            // assertEquals only works with doubles.
            assertEquals(S1Angle(x, boundary_x).degrees(),
                    loop.GetDistanceToBoundary(x).degrees(), kMaxError.degrees())
        }
        if (loop.contains(x)) {
            assertEquals(S1Angle::Zero(), loop.GetDistance(x))
            assertEquals(x, loop.Project(x))
        } else {
            assertEquals(loop.GetDistanceToBoundary(x), loop.GetDistance(x))
            assertEquals(loop.ProjectToBoundary(x), loop.Project(x))
        }
    }

    fun testDistanceMethods() {
        // S2ClosestEdgeQuery is already tested, so just do a bit of sanity checking.

        // The empty and full loops don't have boundaries.
        TestDistanceMethods(*empty, S2Point(0, 1, 0), S2Point())
        TestDistanceMethods(*full, S2Point(0, 1, 0), S2Point())

        // A CCW square around the S2LatLng point (0,0).  Note that because lines of
        // latitude are curved on the sphere, it is not straightforward to project
        // points onto any edge except along the equator.  (The equator is the only
        // line of latitude that is also a geodesic.)
        unique_ptr<S2Loop> square (s2textformat::MakeLoop("-1:-1, -1:1, 1:1, 1:-1"))
        assertTrue(square.isNormalized())

        // A vertex.
        TestDistanceMethods(*square, S2LatLng.fromDegrees(1, -1).toPoint(),
                S2Point())
        // A point on one of the edges.
        TestDistanceMethods(*square, S2LatLng.fromDegrees(0.5, 1).toPoint(),
                S2Point())
        // A point inside the square.
        TestDistanceMethods(*square, S2LatLng.fromDegrees(0, 0.5).toPoint(),
                S2LatLng.fromDegrees(0, 1).toPoint())
        // A point outside the square that projects onto an edge.
        TestDistanceMethods(*square, S2LatLng.fromDegrees(0, -2).toPoint(),
                S2LatLng.fromDegrees(0, -1).toPoint())
        // A point outside the square that projects onto a vertex.
        TestDistanceMethods(*square, S2LatLng.fromDegrees(3, 4).toPoint(),
                S2LatLng.fromDegrees(1, 1).toPoint())
    }

    fun testMakeRegularLoop() {
        S2Point center = S2LatLng . fromDegrees (80, 135).toPoint()
        S1Angle radius = S1Angle ::Degrees(20)
        unique_ptr<S2Loop> loop (S2Loop::MakeRegularLoop(center, radius, 4))

        ASSERT_EQ(4, loop.num_vertices())
        S2Point p0 = loop . vertex (0)
        S2Point p1 = loop . vertex (1)
        S2Point p2 = loop . vertex (2)
        S2Point p3 = loop . vertex (3)
        // Make sure that the radius is correct.
        assertEquals(20.0, S2LatLng(center).GetDistance(S2LatLng(p0)).degrees())
        assertEquals(20.0, S2LatLng(center).GetDistance(S2LatLng(p1)).degrees())
        assertEquals(20.0, S2LatLng(center).GetDistance(S2LatLng(p2)).degrees())
        assertEquals(20.0, S2LatLng(center).GetDistance(S2LatLng(p3)).degrees())
        // Make sure that all angles of the polygon are the same.
        assertEquals(M_PI_2, (p1 - p0).Angle(p3 - p0))
        assertEquals(M_PI_2, (p2 - p1).Angle(p0 - p1))
        assertEquals(M_PI_2, (p3 - p2).Angle(p1 - p2))
        assertEquals(M_PI_2, (p0 - p3).Angle(p2 - p3))
        // Make sure that all edges of the polygon have the same length.
        assertEquals(27.990890717782829,
                S2LatLng(p0).GetDistance(S2LatLng(p1)).degrees())
        assertEquals(27.990890717782829,
                S2LatLng(p1).GetDistance(S2LatLng(p2)).degrees())
        assertEquals(27.990890717782829,
                S2LatLng(p2).GetDistance(S2LatLng(p3)).degrees())
        assertEquals(27.990890717782829,
                S2LatLng(p3).GetDistance(S2LatLng(p0)).degrees())

        // Check actual coordinates. This may change if we switch the algorithm
        // intentionally.
        assertEquals(62.162880741097204, S2LatLng(p0).lat().degrees())
        assertEquals(103.11051028343407, S2LatLng(p0).lng().degrees())
        assertEquals(61.955157772928345, S2LatLng(p1).lat().degrees())
        assertEquals(165.25681963683536, S2LatLng(p1).lng().degrees())
        assertEquals(75.139812547718478, S2LatLng(p2).lat().degrees())
        assertEquals(-119.13042521187423, S2LatLng(p2).lng().degrees())
        assertEquals(75.524190079054392, S2LatLng(p3).lat().degrees())
        assertEquals(26.392175948257943, S2LatLng(p3).lng().degrees())
    }

    TEST(S2LoopShape, Basic)
    {
        unique_ptr<S2Loop> loop = s2textformat ::MakeLoop("0:0, 0:1, 1:0")
        S2Loop::Shape shape (loop.get())
        assertEquals(loop.get(), shape.loop())
        assertEquals(3, shape.num_edges())
        assertEquals(1, shape.num_chains())
        assertEquals(0, shape.chain(0).start)
        assertEquals(3, shape.chain(0).length)
        auto edge2 = shape . edge (2)
        assertEquals("1:0", s2textformat::ToString(edge2.v0))
        assertEquals("0:0", s2textformat::ToString(edge2.v1))
        assertEquals(2, shape.dimension())
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.GetReferencePoint().contained)
    }

    TEST(S2LoopShape, EmptyLoop)
    {
        S2Loop loop (S2Loop.kEmpty)
        S2Loop::Shape shape (&loop)
        assertEquals(0, shape.num_edges())
        assertEquals(0, shape.num_chains())
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.GetReferencePoint().contained)
    }

    TEST(S2LoopShape, FullLoop)
    {
        S2Loop loop (S2Loop.kFull)
        S2Loop::Shape shape (&loop)
        assertEquals(0, shape.num_edges())
        assertEquals(1, shape.num_chains())
        assertFalse(shape.isEmpty())
        assertTrue(shape.isFull())
        assertTrue(shape.GetReferencePoint().contained)
    }

    TEST(S2LoopOwningShape, Ownership)
    {
        // Debug mode builds will catch any memory leak below.
        auto loop = make_unique < S2Loop >(S2Loop.kEmpty)
        S2Loop::OwningShape shape (std::move(loop))
    }
    
 */
}

