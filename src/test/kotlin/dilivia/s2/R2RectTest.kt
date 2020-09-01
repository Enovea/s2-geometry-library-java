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
// Copyright 2012 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

// Author: ericv@google.com (Eric Veach)
//
// Most of the R2Rect methods have trivial implementations in terms of the
// R1Interval class, so most of the testing is done in that unit test.
package dilivia.s2

import dilivia.s2.math.R2Point
import junit.framework.TestCase


@Strictfp
class R2RectTest : TestCase() {

    fun testEmptyRectangles() {
        // Test basic properties of empty rectangles.
        val empty = R2Rect.empty()
        assertTrue(empty.isValid)
        assertTrue(empty.isEmpty)
    }

    fun testConstructorsAndAccessors() {
        // Check various constructors and accessor methods.
        val r = R2Rect(R2Point(0.1, 0.0), R2Point(0.25, 1.0))
        assertEquals(0.1, r.x.lo)
        assertEquals(0.25, r.x.hi)
        assertEquals(0.0, r.y.lo)
        assertEquals(1.0, r.y.hi)

        assertEquals(0.1, r[0][0])
        assertEquals(0.25, r[0][1])
        assertEquals(0.0, r[1][0])
        assertEquals(1.0, r[1][1])

        assertEquals(R1Interval(0.1, 0.25), r.x)
        assertEquals(R1Interval(0, 1), r.y)

        assertEquals(R1Interval(0.1, 0.25), r[0])
        assertEquals(R1Interval(0, 1), r[1])

        val r2 = R2Rect()
        assertTrue(r2.isEmpty)
    }

    fun testFromCenterSize() {
        // FromCenterSize()
        assertTrue(R2Rect.fromCenterSize(R2Point(0.3, 0.5), R2Point(0.2, 0.4)).approxEquals(R2Rect(R2Point(0.2, 0.3), R2Point(0.4, 0.7))))
        assertTrue(R2Rect.fromCenterSize(R2Point(1.0, 0.1), R2Point(0, 2)).approxEquals(R2Rect(R2Point(1.0, -0.9), R2Point(1.0, 1.1))))
    }

    fun testFromPoint() {
        // FromPoint(), FromPointPair()
        val d1 = R2Rect(R2Point(0.1, 0.0), R2Point(0.25, 1.0))
        assertEquals(R2Rect(d1.lo, d1.lo), R2Rect.fromPoint(d1.lo))
        assertEquals(R2Rect(R2Point(0.15, 0.3), R2Point(0.35, 0.9)), R2Rect.fromPointPair(R2Point(0.15, 0.9), R2Point(0.35, 0.3)))
        assertEquals(R2Rect(R2Point(0.12, 0.0), R2Point(0.83, 0.5)), R2Rect.fromPointPair(R2Point(0.83, 0.0), R2Point(0.12, 0.5)))
    }

    fun testSimplePredicates() {
        // GetCenter(), GetVertex(), Contains(R2Point), InteriorContains(R2Point).
        val sw1 = R2Point (0.0, 0.25)
        val ne1 = R2Point (0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)

        assertEquals(R2Point(0.25, 0.5), r1.center)
        assertEquals(R2Point(0.0, 0.25), r1.getVertex(0))
        assertEquals(R2Point(0.5, 0.25), r1.getVertex(1))
        assertEquals(R2Point(0.5, 0.75), r1.getVertex(2))
        assertEquals(R2Point(0.0, 0.75), r1.getVertex(3))
        assertTrue(r1.contains(R2Point(0.2, 0.4)))
        assertFalse(r1.contains(R2Point(0.2, 0.8)))
        assertFalse(r1.contains(R2Point(-0.1, 0.4)))
        assertFalse(r1.contains(R2Point(0.6, 0.1)))
        assertTrue(r1.contains(sw1))
        assertTrue(r1.contains(ne1))
        assertFalse(r1.interiorContains(sw1))
        assertFalse(r1.interiorContains(ne1))

        // Make sure that GetVertex() returns vertices in CCW order.
        for (k in 0 until 4) {
            val a = r1.getVertex(k - 1)
            val b = r1.getVertex(k)
            val c = r1.getVertex(k + 1)
            assertTrue((b - a).ortho().dotProd(c - a) > 0)
        }
    }

    fun testIntervalOperations() {
        // Contains(R2Rect), InteriorContains(R2Rect),
        // Intersects(), InteriorIntersects(), Union(), Intersection().
        //
        // Much more testing of these methods is done in s1interval_test
        // and r1interval_test.

        val empty = R2Rect.empty()
        val sw1 = R2Point (0.0, 0.25)
        val ne1 = R2Point (0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)
        val r1_mid = R2Rect (R2Point(0.25, 0.5), R2Point(0.25, 0.5))
        val r_sw1 = R2Rect(sw1, sw1)
        val r_ne1 = R2Rect(ne1, ne1)

        testIntervalOps(r1, r1_mid, "TTTT", r1, r1_mid)
        testIntervalOps(r1, r_sw1, "TFTF", r1, r_sw1);
        testIntervalOps(r1, r_ne1, "TFTF", r1, r_ne1);

        assertEquals(R2Rect(R2Point(0.0, 0.25), R2Point(0.5, 0.75)), r1)
        testIntervalOps(
                r1,
                R2Rect(R2Point(0.45, 0.1), R2Point(0.75, 0.3)),
                "FFTT",
                R2Rect(R2Point(0.0, 0.1), R2Point(0.75, 0.75)),
                R2Rect(R2Point(0.45, 0.25), R2Point(0.5, 0.3)))
        testIntervalOps(
                r1,
                R2Rect(R2Point(0.5, 0.1), R2Point(0.7, 0.3)),
                "FFTF",
                R2Rect(R2Point(0.0, 0.1), R2Point(0.7, 0.75)),
                R2Rect(R2Point(0.5, 0.25), R2Point(0.5, 0.3)))
        testIntervalOps(r1,
                R2Rect(R2Point(0.45, 0.1), R2Point(0.7, 0.25)),
                "FFTF",
                R2Rect(R2Point(0.0, 0.1), R2Point(0.7, 0.75)),
                R2Rect(R2Point(0.45, 0.25), R2Point(0.5, 0.25)))

        testIntervalOps(
                R2Rect(R2Point(0.1, 0.2), R2Point(0.1, 0.3)),
                R2Rect(R2Point(0.15, 0.7), R2Point(0.2, 0.8)),
                "FFFF",
                R2Rect(R2Point(0.1, 0.2), R2Point(0.2, 0.8)),
                empty)

        // Check that the intersection of two rectangles that overlap in x but not y
        // is valid, and vice versa.
        testIntervalOps(
                R2Rect(R2Point(0.1, 0.2), R2Point(0.4, 0.5)),
                R2Rect(R2Point(0, 0), R2Point(0.2, 0.1)),
                "FFFF",
                R2Rect(R2Point(0, 0), R2Point(0.4, 0.5)), empty)
        testIntervalOps(
                R2Rect(R2Point(0, 0), R2Point(0.1, 0.3)),
                R2Rect(R2Point(0.2, 0.1), R2Point(0.3, 0.4)),
                "FFFF",
                R2Rect(R2Point(0, 0), R2Point(0.3, 0.4)), empty)
    }

    fun testAddPoint() {
        // AddPoint()
        val sw1 = R2Point (0.0, 0.25)
        val ne1 = R2Point (0.5, 0.75)
        val r1 = R2Rect(sw1, ne1)

        var r2 = R2Rect.empty()
        r2 = r2.addPoint(R2Point(0.0, 0.25))
        r2 = r2.addPoint(R2Point(0.5, 0.25))
        r2 = r2.addPoint(R2Point(0.0, 0.75))
        r2 = r2.addPoint(R2Point(0.1, 0.4))
        assertEquals(r1, r2);
    }

    fun testProject() {
        val r1 = R2Rect(R1Interval(0.0, 0.5), R1Interval(0.25, 0.75))

        assertEquals(R2Point(0.0, 0.25), r1.project(R2Point(-0.01, 0.24)))
        assertEquals(R2Point(0.0, 0.48), r1.project(R2Point(-5.0, 0.48)))
        assertEquals(R2Point(0.0, 0.75), r1.project(R2Point(-5.0, 2.48)))
        assertEquals(R2Point(0.19, 0.75), r1.project(R2Point(0.19, 2.48)))
        assertEquals(R2Point(0.5, 0.75), r1.project(R2Point(6.19, 2.48)))
        assertEquals(R2Point(0.5, 0.53), r1.project(R2Point(6.19, 0.53)))
        assertEquals(R2Point(0.5, 0.25), r1.project(R2Point(6.19, -2.53)))
        assertEquals(R2Point(0.33, 0.25), r1.project(R2Point(0.33, -2.53)))
        assertEquals(R2Point(0.33, 0.37), r1.project(R2Point(0.33, 0.37)))
    }

    fun testExpanded() {
        // Expanded()
        assertTrue(R2Rect.empty().expanded(R2Point(0.1, 0.3)).isEmpty);
        assertTrue(R2Rect.empty().expanded(R2Point(-0.1, -0.3)).isEmpty);
        assertTrue(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, 0.3)).approxEquals(R2Rect(R2Point(0.1, 0.1), R2Point(0.4, 1.0))));
        assertTrue(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(-0.1, 0.3)).isEmpty);
        assertTrue(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, -0.2)).isEmpty);
        assertTrue(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(R2Point(0.1, -0.1)).approxEquals(R2Rect(R2Point(0.1, 0.5), R2Point(0.4, 0.6))));
        assertTrue(R2Rect(R2Point(0.2, 0.4), R2Point(0.3, 0.7)).expanded(0.1).approxEquals(R2Rect(R2Point(0.1, 0.3), R2Point(0.4, 0.8))));
    }

    companion object {
        fun testIntervalOps(x: R2Rect, y: R2Rect, expectedRexion: String, expectedUnion: R2Rect, expectedIntersection: R2Rect) {
            // Test all of the interval operations on the given pair of intervals.
            // "expected_rexion" is a sequence of "T" and "F" characters corresponding
            // to the expected results of Contains(), InteriorContains(), Intersects(),
            // and InteriorIntersects() respectively.

            assertEquals(expectedRexion[0] == 'T', x.contains(y));
            assertEquals(expectedRexion[1] == 'T', x.interiorContains(y));
            assertEquals(expectedRexion[2] == 'T', x.intersects(y));
            assertEquals(expectedRexion[3] == 'T', x.interiorIntersects(y));

            assertEquals(x.union(y) == x, x.contains(y));
            assertEquals(!x.intersection(y).isEmpty, x.intersects(y));

            assertEquals(expectedUnion, x.union(y));
            assertEquals(expectedIntersection, x.intersection(y));

            var r = x;
            r = r.addRect(y)
            assertEquals(expectedUnion, r);
            if (y.size == R2Point(0, 0)) {
                r = x;
                r = r.addPoint(y.lo)
                assertEquals(expectedUnion, r)
            }
        }
    }
}