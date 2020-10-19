/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package dilivia.s2.region

import dilivia.s2.S2.M_PI
import dilivia.s2.S2.M_PI_2
import dilivia.s2.*
import dilivia.s2.region.S2Loop.S2PointLoopSpan
import dilivia.s2.region.S2LoopMeasures.LoopOrder
import kotlin.math.*

class S2LoopMeasuresTest : S2GeometryTestCase() {

    // Some standard loops to use in the tests (see descriptions below).
    // The full loop is represented as a loop with no vertices.
    val full = S2PointLoopSpan()

    // A degenerate loop in the shape of a "V".
    val v_loop = S2PointLoopSpan(parseVertices("5:1, 0:2, 5:3, 0:2"))

    // The northern hemisphere, defined using two pairs of antipodal points.
    val north_hemi = S2PointLoopSpan(parseVertices("0:-180, 0:-90, 0:0, 0:90"))

    // The northern hemisphere, defined using three points 120 degrees apart.
    val north_hemi3 = S2PointLoopSpan(parseVertices("0:-180, 0:-60, 0:60"))

    // The western hemisphere, defined using two pairs of antipodal points.
    val west_hemi = S2PointLoopSpan(parseVertices("0:-180, -90:0, 0:0, 90:0"))

    // The eastern hemisphere, defined using two pairs of antipodal points.
    val east_hemi = S2PointLoopSpan(parseVertices("90:0, 0:0, -90:0, 0:-180"))

    // A spiral stripe that slightly over-wraps the equator.
    val candy_cane = S2PointLoopSpan(parseVertices("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70"))

    // A completely degenerate triangle along the equator that Sign()
    // considers to be CCW.
    val line_triangle = S2PointLoopSpan(parseVertices("0:1, 0:2, 0:3"))

    // A nearly-degenerate CCW chevron near the equator with very long sides
    // (about 80 degrees).  Its area is less than 1e-640, which is too small
    // to represent in double precision.
    val skinny_chevron = S2PointLoopSpan(parseVertices("0:0, -1e-320:80, 0:1e-320, 1e-320:80"))

    // A loop where the same vertex appears three times.
    val three_leaf_clover = S2PointLoopSpan(parseVertices("0:0, -3:3, 3:3, 0:0, 3:0, 3:-3, 0:0, -3:-3, -3:0"))

    // A loop with groups of 3 or more vertices in a straight line.
    val tessellated_loop = S2PointLoopSpan(parseVertices("10:34, 5:34, 0:34, -10:34, -10:36, -5:36, 0:36, 10:36"))

    // Given a string where each character "ch" represents a vertex (such as
// "abac"), returns a vector of S2Points of the form (ch, 0, 0).  Note that
// these points are not unit length and therefore are not suitable for general
// use, however they are useful for testing certain functions below.
    fun makeTestLoop(loop_str: String): List<S2Point> {
        val loop = mutableListOf<S2Point>();
        for (ch in loop_str) {
            loop.add(S2Point(ch.toInt(), 0, 0));
        }
        return loop;
    }

    // Given a loop whose vertices are represented as characters (such as "abcd" or
// "abccb"), verify that S2::PruneDegeneracies() yields the loop "expected".
    fun testPruneDegeneracies(input_str: String, expected_str: String) {
        val input = makeTestLoop(input_str)
        val pruned = S2LoopMeasures.pruneDegeneracies(S2PointLoopSpan(input))
        var actual_str: String = ""
        for (p in pruned) {
            actual_str += p[0].toChar()
        }
        assertEquals(expected_str, actual_str)
    }

    fun testPruneDegeneraciesAllDegeneracies() {
        testPruneDegeneracies("", "");
        testPruneDegeneracies("a", "");
        testPruneDegeneracies("aaaaa", "");
        testPruneDegeneracies("ab", "");
        testPruneDegeneracies("abb", "");
        testPruneDegeneracies("aab", "");
        testPruneDegeneracies("aba", "");
        testPruneDegeneracies("abba", "");
        testPruneDegeneracies("abcb", "");
        testPruneDegeneracies("abcba", "");
        testPruneDegeneracies("abcdcdedefedcbcdcb", "");
    }

    fun testPruneDegeneraciesSomeDegeneracies() {
        testPruneDegeneracies("abc", "abc");
        testPruneDegeneracies("abca", "abc");
        testPruneDegeneracies("abcc", "abc");
        testPruneDegeneracies("abccaa", "abc");
        testPruneDegeneracies("aabbcc", "abc");
        testPruneDegeneracies("abcdedca", "abc");
        testPruneDegeneracies("abcbabcbcdc", "abc");
        testPruneDegeneracies("xyzabcazy", "abc");
        testPruneDegeneracies("xxyyzzaabbccaazzyyxx", "abc");
    }

    // Given a loop whose vertices are represented as characters (such as "abcd" or
// "abccb"), verify that S2::GetCanonicalLoopOrder returns the given result.
    fun testCanonicalLoopOrder(input_str: String, expected_order: LoopOrder) {
        assertEquals(expected_order, S2LoopMeasures.getCanonicalLoopOrder(S2PointLoopSpan(makeTestLoop(input_str))))
    }

    fun testGetCanonicalLoopOrderAllDegeneracies() {
        testCanonicalLoopOrder("", LoopOrder(0, 1))
        testCanonicalLoopOrder("a", LoopOrder(0, 1))
        testCanonicalLoopOrder("aaaaa", LoopOrder(0, 1))
        testCanonicalLoopOrder("ba", LoopOrder(1, 1))
        testCanonicalLoopOrder("bab", LoopOrder(1, 1))
        testCanonicalLoopOrder("cbab", LoopOrder(2, 1))
        testCanonicalLoopOrder("bacbcab", LoopOrder(8, -1))
    }

    fun testGetPerimeterEmpty() {
        assertEquals(S1Angle.zero, S2LoopMeasures.getPerimeter(S2PointLoopSpan()))
    }

    fun testGetPerimeterOctant() {
        val loop = parseVertices("0:0, 0:90, 90:0");
        assertEquals(3 * M_PI_2, S2LoopMeasures.getPerimeter(S2PointLoopSpan(loop)).radians)
    }

    fun testGetPerimeterMoreThanTwoPi() {
        // Make sure that GetPerimeter doesn't use S1ChordAngle, which can only
        // represent distances up to 2*Pi.
        val loop = parseVertices("0:0, 0:90, 0:180, 90:0, 0:-90");
        assertEquals(5 * M_PI_2, S2LoopMeasures.getPerimeter(S2PointLoopSpan(loop)).radians)
    }


    fun testAreaConsistentWithCurvature(loop: S2PointLoopSpan) {
        // Check that the area computed using GetArea() is consistent with the loop
        // curvature.  According to the Gauss-Bonnet theorem, the area of the loop
        // equals 2*Pi minus its curvature.
        val area = S2LoopMeasures.getArea(loop)
        val gauss_area = 2 * M_PI - S2LoopMeasures.getCurvature(loop)
        // The error bound below is sufficient for current tests but not guaranteed.
        assertTrue("Failed loop: \n$loop\nArea = $area, Gauss Area = $gauss_area", abs(area - gauss_area) <= 1e-14)

    }

    fun testGetAreaConsistentWithCurvature() {
        testAreaConsistentWithCurvature(full);
        testAreaConsistentWithCurvature(north_hemi);
        testAreaConsistentWithCurvature(north_hemi3);
        testAreaConsistentWithCurvature(west_hemi);
        testAreaConsistentWithCurvature(east_hemi);
        testAreaConsistentWithCurvature(candy_cane);
        testAreaConsistentWithCurvature(line_triangle);
        testAreaConsistentWithCurvature(skinny_chevron);
        testAreaConsistentWithCurvature(three_leaf_clover);
        testAreaConsistentWithCurvature(tessellated_loop);
    }

    private val kMaxVertices = 6

    fun testGetAreaConsistentWithOrientation() {
        // Test that GetArea() returns an area near 0 for degenerate loops that
        // contain almost no points, and an area near 4*Pi for degenerate loops that
        // contain almost all points.

        repeat(50) { i ->
            val num_vertices = 3 + S2Random.randomInt(kMaxVertices - 3 + 1)
            // Repeatedly choose N vertices that are exactly on the equator until we
            // find some that form a valid loop.
            val loop = mutableListOf<S2Point>()
            do {
                loop.clear();
                repeat(num_vertices) {
                    // We limit longitude to the range [0, 90] to ensure that the loop is
                    // degenerate (as opposed to following the entire equator).
                    loop.add(S2LatLng.fromRadians(0.0, S2Random.randomDouble() * M_PI_2).toPoint())
                }
            } while (!S2Loop(loop, check = false).isValid())
            val ccw = S2LoopMeasures.isNormalized(S2PointLoopSpan(loop))
            // The error bound is sufficient for current tests but not guaranteed.
            assertEquals("Failed loop $i: $loop", if (ccw) 0.0 else 4.0 * M_PI, S2LoopMeasures.getArea(S2PointLoopSpan(loop)), 1e-14)

            assertEquals(!ccw, S2Loop(loop).contains(S2Point(0, 0, 1)));
        }
    }

    fun testGetAreaAccuracy() {
        // TODO(ericv): Test that GetArea() has an accuracy significantly better
        // than 1e-15 on loops whose area is small.
    }

    val kMaxDist = 1e-6
    fun testGetAreaAndCentroid() {
        assertEquals(4 * M_PI, S2LoopMeasures.getArea(full));
        assertEquals(S2Point(0, 0, 0), S2LoopMeasures.getCentroid(full))

        assertEquals(S2LoopMeasures.getArea(north_hemi), 2 * M_PI);
        assertEquals(2 * M_PI, S2LoopMeasures.getArea(east_hemi), 1e-12);

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices.  Then check that the area and centroid are
        // correct.
        repeat(50) { iter ->
            // Choose a coordinate frame for the spherical cap.
          val (x, y, z) = S2Random.randomFrame()

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            val height = 2 * S2Random.randomDouble()
            val phi = asin(1 - height)
            var max_dtheta = 2 * acos(tan(abs(phi)) / tan(abs(phi) + kMaxDist))
            max_dtheta = min(M_PI, max_dtheta);  // At least 3 vertices.

            val loop = mutableListOf<S2Point>()
            var theta = 0.0
            while (theta < 2 * M_PI) {
                loop.add(cos(theta) * cos(phi) * x + sin(theta) * cos(phi) * y + sin(phi) * z)
                theta += S2Random.randomDouble() * max_dtheta
            }
            val area = S2LoopMeasures.getArea(S2PointLoopSpan(loop))
            val centroid = S2LoopMeasures.getCentroid(S2PointLoopSpan(loop))
            val expected_area = 2 * M_PI * height
            assertTrue(abs(area - expected_area) <= 2 * M_PI * kMaxDist)
            val expected_centroid = expected_area * (1 - 0.5 * height) * z
            assertTrue((centroid - expected_centroid).norm() <= 2 * kMaxDist)
        }
    }

    fun expectSameOrder(loop1: S2PointLoopSpan, order1: LoopOrder, loop2: S2PointLoopSpan, order2: LoopOrder) {
        assertEquals(loop1.size, loop2.size)
        var i1 = order1.first
        var i2 = order2.first
        val dir1 = order1.dir
        val dir2 = order2.dir
        var n = loop1.size
        while (--n >= 0) {
            assertEquals("$order1 vs. $order2", loop1[i1], loop2[i2])
            i1 += dir1;
            i2 += dir2;
        }
    }

    // Check that the curvature is *identical* when the vertex order is
// rotated, and that the sign is inverted when the vertices are reversed.
    fun checkCurvatureInvariants(loop_in: S2PointLoopSpan) {
        val order_in = S2LoopMeasures.getCanonicalLoopOrder(loop_in)
        var loop: S2PointLoopSpan = S2PointLoopSpan(loop_in.points)
        val expected = S2LoopMeasures.getCurvature(loop_in);
        for (i in 0 until loop_in.size) {
            loop = S2PointLoopSpan(loop.points.reversed())
            assertEquals(if (expected == 2 * M_PI) expected else -expected, S2LoopMeasures.getCurvature(loop))
            expectSameOrder(loop_in, order_in, loop, S2LoopMeasures.getCanonicalLoopOrder(loop))
            loop = S2PointLoopSpan(loop.points.reversed())
            loop = loop.rotated(1)
            assertEquals(expected, S2LoopMeasures.getCurvature(loop));
            expectSameOrder(loop_in, order_in, loop, S2LoopMeasures.getCanonicalLoopOrder(loop));
        }
    }

    fun testGetCurvature() {
        assertEquals(-2 * M_PI, S2LoopMeasures.getCurvature(full))

        assertEquals(2 * M_PI, S2LoopMeasures.getCurvature(v_loop))
        checkCurvatureInvariants(v_loop);

        // This curvature should be computed exactly.
        assertEquals(0.0, S2LoopMeasures.getCurvature(north_hemi3))
        checkCurvatureInvariants(north_hemi3)

        assertEquals(0.0, S2LoopMeasures.getCurvature(west_hemi), 1e-15)
        checkCurvatureInvariants(west_hemi)

        // We don't have an easy way to estimate the curvature of these loops, but
        // we can still check that the expected invariants hold.
        checkCurvatureInvariants(candy_cane)
        checkCurvatureInvariants(three_leaf_clover)

        assertEquals(2 * M_PI, S2LoopMeasures.getCurvature(line_triangle), 1e-15)
        checkCurvatureInvariants(line_triangle)

        assertEquals(2 * M_PI, S2LoopMeasures.getCurvature(skinny_chevron), 1e-15)
        checkCurvatureInvariants(skinny_chevron)

        // Build a narrow spiral loop starting at the north pole.  This is designed
        // to test that the error in GetCurvature is linear in the number of
        // vertices even when the partial sum of the curvatures gets very large.
        // The spiral consists of two "arms" defining opposite sides of the loop.
        // This is a pathological loop that contains many long parallel edges.
        val kArmPoints = 10000;    // Number of vertices in each "arm"
        val kArmRadius = 0.01;  // Radius of spiral.
        val spiral = mutableListOf<S2Point>()
        repeat(2 * kArmPoints) { spiral.add(S2Point(0, 0, 1)) }
        for (i in 0 until kArmPoints) {
            val angle = (2 * M_PI / 3) * i;
            val x = cos(angle);
            val y = sin(angle);
            val r1 = i * kArmRadius / kArmPoints;
            val r2 = (i + 1.5) * kArmRadius / kArmPoints;
            spiral[kArmPoints - i - 1] = S2Point(r1 * x, r1 * y, 1.0).normalize();
            spiral[kArmPoints + i] = S2Point(r2 * x, r2 * y, 1.0).normalize();
        }

        // Check that GetCurvature() is consistent with GetArea() to within the
        // error bound of the former.  We actually use a tiny fraction of the
        // worst-case error bound, since the worst case only happens when all the
        // roundoff errors happen in the same direction and this test is not
        // designed to achieve that.  The error in GetArea() can be ignored for the
        // purposes of this test since it is generally much smaller.
        assertEquals(2 * M_PI - S2LoopMeasures.getArea(S2PointLoopSpan(spiral)), S2LoopMeasures.getCurvature(S2PointLoopSpan(spiral)),
                0.01 * S2LoopMeasures.getCurvatureMaxError(S2PointLoopSpan(spiral)))
    }

}
