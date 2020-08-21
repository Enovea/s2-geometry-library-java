/**
 * This project is a kotlin port of the Google s2 geometry library: https://github.com/google/s2geometry.git
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
package dilivia.s2

import com.google.common.geometry.GeometryTestCase
import com.google.common.geometry.S2.*
import kotlin.math.abs

@Strictfp
class S1IntervalTest : GeometryTestCase() {

    // "Quadrants" are numbered as follows:
    // quad1 == [0, Pi/2]
    // quad2 == [Pi/2, Pi]
    // quad3 == [-Pi, -Pi/2]
    // quad4 == [-Pi/2, 0]


    val empty = S1Interval.empty
    val full = S1Interval.full

    val zero = S1Interval(0.0, 0.0)
    val pi2 = S1Interval(M_PI_2, M_PI_2)
    val pi = S1Interval(M_PI, M_PI)
    val mipi = S1Interval(-M_PI, -M_PI)
    val mipi2 = S1Interval(-M_PI_2, -M_PI_2)

    val quad1 = S1Interval(0.0, M_PI_2)
    val quad2 = S1Interval(M_PI_2, -M_PI)
    val quad3 = S1Interval(M_PI, -M_PI_2)
    val quad4 = S1Interval(-M_PI_2, 0.0)

    val quad12 = S1Interval(0.0, -M_PI)
    val quad23 = S1Interval(M_PI_2, -M_PI_2) // inverted
    val quad34 = S1Interval(-M_PI, 0.0)
    val quad41 = S1Interval(-M_PI_2, M_PI_2)

    val quad123 = S1Interval(0.0, -M_PI_2)
    val quad234 = S1Interval(M_PI_2, 0.0)
    val quad341 = S1Interval(M_PI, M_PI_2)
    val quad412 = S1Interval(-M_PI_2, -M_PI)

    val mid12 = S1Interval(M_PI_2 - 0.01, M_PI_2 + 0.02)
    val mid23 = S1Interval(M_PI - 0.01, -M_PI + 0.02)
    val mid34 = S1Interval(-M_PI_2 - 0.01, -M_PI_2 + 0.02)
    val mid41 = S1Interval(-0.01, 0.02)

    val quad1lo = S1Interval(quad12.lo, mid41.hi)
    val quad2lo = S1Interval(quad23.lo, mid12.hi)
    val quad2hi = S1Interval(mid23.lo, quad12.hi)
    val quad3hi = S1Interval(mid34.lo, quad23.hi)
    val quad12eps = S1Interval(quad12.lo, mid23.hi)
    val quadeps12 = S1Interval(mid41.lo, quad12.hi)
    val quad123eps = S1Interval(quad12.lo, mid34.hi)
    val quadeps123 = S1Interval(mid41.lo, quad23.hi)
    val quad23eps = S1Interval(quad23.lo, mid34.hi)
    val quadeps23 = S1Interval(mid12.lo, quad23.hi)
    val quad412eps = S1Interval(mid34.lo, quad12.hi)

    fun testConstructorsAndAccessors() {
        // Spot-check the constructors and accessors.
        assertEquals(quad12.lo, 0.0)
        assertEquals(quad12.hi, M_PI)
        assertEquals(quad34[0], M_PI)
        assertEquals(quad34[1], 0.0)
        assertEquals(quad34.bounds, R2Vector(M_PI, 0.0))
        assertEquals(pi.lo, M_PI)
        assertEquals(pi.hi, M_PI)

        // Check that [-Pi, -Pi] is normalized to [Pi, Pi].
        assertEquals(mipi.lo, M_PI)
        assertEquals(mipi.hi, M_PI)
        assertEquals(quad23.lo, M_PI_2)
        assertEquals(quad23.hi, -M_PI_2)

        // Check that the default S1Interval is identical to Empty().
        val defaultEmpty = S1Interval()
        assertTrue(defaultEmpty.isValid)
        assertTrue(defaultEmpty.isEmpty)
        assertEquals(empty.lo, defaultEmpty.lo)
        assertEquals(empty.hi, defaultEmpty.hi)
    }

    fun testSimplePredicates() {
        // is_valid(), is_empty(), is_full(), is_inverted()
        assertTrue(zero.isValid && !zero.isEmpty && !zero.isFull)
        assertTrue(empty.isValid && empty.isEmpty && !empty.isFull)
        assertTrue(empty.isInverted)
        assertTrue(full.isValid && !full.isEmpty && full.isFull)
        assertTrue(!quad12.isEmpty && !quad12.isFull && !quad12.isInverted)
        assertTrue(!quad23.isEmpty && !quad23.isFull && quad23.isInverted)
        assertTrue(pi.isValid && !pi.isEmpty && !pi.isInverted)
        assertTrue(mipi.isValid && !mipi.isEmpty && !mipi.isInverted)
    }

    fun testAlmostEmptyOrFull() {
        // Test that rounding errors don't cause intervals that are almost empty or
        // full to be considered empty or full.  The following value is the greatest
        // representable value less than Pi.
        val kAlmostPi = M_PI - 2 * DBL_EPSILON
        assertFalse(S1Interval(-kAlmostPi, M_PI).isFull)
        assertFalse(S1Interval(-M_PI, kAlmostPi).isFull)
        assertFalse(S1Interval(M_PI, -kAlmostPi).isEmpty)
        assertFalse(S1Interval(kAlmostPi, -M_PI).isEmpty)
    }

    fun testGetCenter() {
        assertEquals(quad12.center, M_PI_2)
        assertEquals(S1Interval(3.1, 2.9).center, 3.0 - M_PI)
        assertEquals(S1Interval(-2.9, -3.1).center, M_PI - 3.0)
        assertEquals(S1Interval(2.1, -2.1).center, M_PI)
        assertEquals(pi.center, M_PI)
        assertEquals(mipi.center, M_PI)
        assertEquals(abs(quad23.center), M_PI)
        assertEquals(quad123.center, 0.75 * M_PI)
    }

    fun testGetLength() {
        assertEquals(quad12.length, M_PI)
        assertEquals(pi.length, 0.0)
        assertEquals(mipi.length, 0.0)
        assertEquals(quad123.length, 1.5 * M_PI)
        assertEquals(abs(quad23.length), M_PI)
        assertEquals(full.length, 2 * M_PI)
        assertTrue(empty.length < 0)
    }

    fun testComplement() {
        assertTrue(empty.complement.isFull)
        assertTrue(full.complement.isEmpty)
        assertTrue(pi.complement.isFull)
        assertTrue(mipi.complement.isFull)
        assertTrue(zero.complement.isFull)
        assertTrue(quad12.complement.approxEquals(quad34))
        assertTrue(quad34.complement.approxEquals(quad12))
        assertTrue(quad123.complement.approxEquals(quad4))
    }

    fun testContains() {
        // Contains(double), InteriorContains(double)
        assertTrue(!empty.contains(0) && !empty.contains(M_PI) && !empty.contains(-M_PI))
        assertTrue(!empty.interiorContains(M_PI) && !empty.interiorContains(-M_PI))
        assertTrue(full.contains(0) && full.contains(M_PI) && full.contains(-M_PI))
        assertTrue(full.interiorContains(M_PI) && full.interiorContains(-M_PI))
        assertTrue(quad12.contains(0) && quad12.contains(M_PI) && quad12.contains(-M_PI))
        assertTrue(quad12.interiorContains(M_PI_2) && !quad12.interiorContains(0))
        assertTrue(!quad12.interiorContains(M_PI) && !quad12.interiorContains(-M_PI))
        assertTrue(quad23.contains(M_PI_2) && quad23.contains(-M_PI_2))
        assertTrue(quad23.contains(M_PI) && quad23.contains(-M_PI))
        assertTrue(!quad23.contains(0))
        assertTrue(!quad23.interiorContains(M_PI_2) && !quad23.interiorContains(-M_PI_2))
        assertTrue(quad23.interiorContains(M_PI) && quad23.interiorContains(-M_PI))
        assertTrue(!quad23.interiorContains(0))
        assertTrue(pi.contains(M_PI) && pi.contains(-M_PI) && !pi.contains(0))
        assertTrue(!pi.interiorContains(M_PI) && !pi.interiorContains(-M_PI))
        assertTrue(mipi.contains(M_PI) && mipi.contains(-M_PI) && !mipi.contains(0))
        assertTrue(!mipi.interiorContains(M_PI) && !mipi.interiorContains(-M_PI))
        assertTrue(zero.contains(0) && !zero.interiorContains(0))
    }

    private fun testIntervalOps(x: S1Interval, y: S1Interval, expectedRelation: String,
                                expectedUnion: S1Interval, expectedIntersection: S1Interval) {
        // Test all of the interval operations on the given pair of intervals.
        // "expected_relation" is a sequence of "T" and "F" characters corresponding
        // to the expected results of Contains(), InteriorContains(), Intersects(),
        // and InteriorIntersects() respectively.
        assertEquals(x.contains(y), expectedRelation[0] == 'T')
        assertEquals(x.interiorContains(y), expectedRelation[1] == 'T')
        assertEquals(x.intersects(y), expectedRelation[2] == 'T')
        assertEquals(x.interiorIntersects(y), expectedRelation[3] == 'T')

        // bounds() returns a const reference to a member variable, so we need to
        // make a copy when invoking it on a temporary object.
        println("-----------")
        println("x = $x")
        println("y = $y")
        println("x U y = ${x.union(y)}")
        assertEquals(expectedUnion, x.union(y))
        assertEquals(expectedIntersection, x.intersection(y))
        assertEquals(x.contains(y), x.union(y) === x)
        assertEquals(x.intersects(y), !x.intersection(y).isEmpty)
        if (y.lo == y.hi) {
            val r = x.addPoint(y.lo)
            assertEquals(expectedUnion, r)
        }
    }


    fun testIntervalOps() {
        // Contains(S1Interval), InteriorContains(S1Interval),
        // Intersects(), InteriorIntersects(), Union(), Intersection()
        testIntervalOps(empty, empty, "TTFF", empty, empty)
        testIntervalOps(empty, full, "FFFF", full, empty)
        testIntervalOps(empty, zero, "FFFF", zero, empty)
        testIntervalOps(empty, pi, "FFFF", pi, empty)
        testIntervalOps(empty, mipi, "FFFF", mipi, empty)

        testIntervalOps(full, empty, "TTFF", full, empty)
        testIntervalOps(full, full, "TTTT", full, full)
        testIntervalOps(full, zero, "TTTT", full, zero)
        testIntervalOps(full, pi, "TTTT", full, pi)
        testIntervalOps(full, mipi, "TTTT", full, mipi)
        testIntervalOps(full, quad12, "TTTT", full, quad12)
        testIntervalOps(full, quad23, "TTTT", full, quad23)

        testIntervalOps(zero, empty, "TTFF", zero, empty)
        testIntervalOps(zero, full, "FFTF", full, zero)
        testIntervalOps(zero, zero, "TFTF", zero, zero)
        testIntervalOps(zero, pi, "FFFF", S1Interval(0.0, M_PI), empty)
        testIntervalOps(zero, pi2, "FFFF", quad1, empty)
        testIntervalOps(zero, mipi, "FFFF", quad12, empty)
        testIntervalOps(zero, mipi2, "FFFF", quad4, empty)
        testIntervalOps(zero, quad12, "FFTF", quad12, zero)
        testIntervalOps(zero, quad23, "FFFF", quad123, empty)

        testIntervalOps(pi2, empty, "TTFF", pi2, empty)
        testIntervalOps(pi2, full, "FFTF", full, pi2)
        testIntervalOps(pi2, zero, "FFFF", quad1, empty)
        testIntervalOps(pi2, pi, "FFFF", S1Interval(M_PI_2, M_PI), empty)
        testIntervalOps(pi2, pi2, "TFTF", pi2, pi2)
        testIntervalOps(pi2, mipi, "FFFF", quad2, empty)
        testIntervalOps(pi2, mipi2, "FFFF", quad23, empty)
        testIntervalOps(pi2, quad12, "FFTF", quad12, pi2)
        testIntervalOps(pi2, quad23, "FFTF", quad23, pi2)

        testIntervalOps(pi, empty, "TTFF", pi, empty)
        testIntervalOps(pi, full, "FFTF", full, pi)
        testIntervalOps(pi, zero, "FFFF", S1Interval(M_PI, 0.0), empty)
        testIntervalOps(pi, pi, "TFTF", pi, pi)
        testIntervalOps(pi, pi2, "FFFF", S1Interval(M_PI_2, M_PI), empty)
        testIntervalOps(pi, mipi, "TFTF", pi, pi)
        testIntervalOps(pi, mipi2, "FFFF", quad3, empty)
        testIntervalOps(pi, quad12, "FFTF", S1Interval(0.0, M_PI), pi)
        testIntervalOps(pi, quad23, "FFTF", quad23, pi)

        testIntervalOps(mipi, empty, "TTFF", mipi, empty)
        testIntervalOps(mipi, full, "FFTF", full, mipi)
        testIntervalOps(mipi, zero, "FFFF", quad34, empty)
        testIntervalOps(mipi, pi, "TFTF", mipi, mipi)
        testIntervalOps(mipi, pi2, "FFFF", quad2, empty)
        testIntervalOps(mipi, mipi, "TFTF", mipi, mipi)
        testIntervalOps(mipi, mipi2, "FFFF", S1Interval(-M_PI, -M_PI_2), empty)
        testIntervalOps(mipi, quad12, "FFTF", quad12, mipi)
        testIntervalOps(mipi, quad23, "FFTF", quad23, mipi)

        testIntervalOps(quad12, empty, "TTFF", quad12, empty)
        testIntervalOps(quad12, full, "FFTT", full, quad12)
        testIntervalOps(quad12, zero, "TFTF", quad12, zero)
        testIntervalOps(quad12, pi, "TFTF", quad12, pi)
        testIntervalOps(quad12, mipi, "TFTF", quad12, mipi)
        testIntervalOps(quad12, quad12, "TFTT", quad12, quad12)
        testIntervalOps(quad12, quad23, "FFTT", quad123, quad2)
        testIntervalOps(quad12, quad34, "FFTF", full, quad12)

        testIntervalOps(quad23, empty, "TTFF", quad23, empty)
        testIntervalOps(quad23, full, "FFTT", full, quad23)
        testIntervalOps(quad23, zero, "FFFF", quad234, empty)
        testIntervalOps(quad23, pi, "TTTT", quad23, pi)
        testIntervalOps(quad23, mipi, "TTTT", quad23, mipi)
        testIntervalOps(quad23, quad12, "FFTT", quad123, quad2)
        testIntervalOps(quad23, quad23, "TFTT", quad23, quad23)
        testIntervalOps(quad23, quad34, "FFTT", quad234, S1Interval(-M_PI, -M_PI_2))

        testIntervalOps(quad1, quad23, "FFTF", quad123, S1Interval(M_PI_2, M_PI_2))
        testIntervalOps(quad2, quad3, "FFTF", quad23, mipi)
        testIntervalOps(quad3, quad2, "FFTF", quad23, pi)
        testIntervalOps(quad2, pi, "TFTF", quad2, pi)
        testIntervalOps(quad2, mipi, "TFTF", quad2, mipi)
        testIntervalOps(quad3, pi, "TFTF", quad3, pi)
        testIntervalOps(quad3, mipi, "TFTF", quad3, mipi)

        testIntervalOps(quad12, mid12, "TTTT", quad12, mid12)
        testIntervalOps(mid12, quad12, "FFTT", quad12, mid12)

        testIntervalOps(quad12, mid23, "FFTT", quad12eps, quad2hi)
        testIntervalOps(mid23, quad12, "FFTT", quad12eps, quad2hi)

        // This test checks that the union of two disjoint intervals is the smallest
        // interval that contains both of them.  Note that the center of "mid34"
        // slightly CCW of -Pi/2 so that there is no ambiguity about the result.
        testIntervalOps(quad12, mid34, "FFFF", quad412eps, empty)
        testIntervalOps(mid34, quad12, "FFFF", quad412eps, empty)

        testIntervalOps(quad12, mid41, "FFTT", quadeps12, quad1lo)
        testIntervalOps(mid41, quad12, "FFTT", quadeps12, quad1lo)

        testIntervalOps(quad23, mid12, "FFTT", quadeps23, quad2lo)
        testIntervalOps(mid12, quad23, "FFTT", quadeps23, quad2lo)
        testIntervalOps(quad23, mid23, "TTTT", quad23, mid23)
        testIntervalOps(mid23, quad23, "FFTT", quad23, mid23)
        testIntervalOps(quad23, mid34, "FFTT", quad23eps, quad3hi)
        testIntervalOps(mid34, quad23, "FFTT", quad23eps, quad3hi)
        testIntervalOps(quad23, mid41, "FFFF", quadeps123, empty)
        testIntervalOps(mid41, quad23, "FFFF", quadeps123, empty)
    }

    fun testAddPoint() {
        assertEquals(zero, empty.addPoint(0))
        assertEquals(pi, empty.addPoint(M_PI))
        assertEquals(mipi, empty.addPoint(-M_PI))
        assertEquals(pi, empty.addPoint(M_PI).addPoint(-M_PI))
        assertEquals(mipi, empty.addPoint(-M_PI).addPoint(M_PI))
        assertEquals(mid12, empty.addPoint(mid12.lo).addPoint(mid12.hi))
        assertEquals(mid23, empty.addPoint(mid23.lo).addPoint(mid23.hi))
        assertEquals(quad123, quad1.addPoint(-0.9 * M_PI).addPoint(-M_PI_2))
        assertTrue(full.addPoint(0).isFull)
        assertTrue(full.addPoint(M_PI).isFull)
        assertTrue(full.addPoint(-M_PI).isFull)
    }

    fun testProject() {
        var r = S1Interval(-M_PI, -M_PI)
        assertEquals(M_PI, r.project(-M_PI))
        assertEquals(M_PI, r.project(0))
        r = S1Interval(0.0, M_PI)
        assertEquals(0.1, r.project(0.1))
        assertEquals(0.0, r.project(-M_PI_2 + 1e-15))
        assertEquals(M_PI, r.project(-M_PI_2 - 1e-15))
        r = S1Interval(M_PI - 0.1, -M_PI + 0.1)
        assertEquals(M_PI, r.project(M_PI))
        assertEquals(M_PI - 0.1, r.project(1e-15))
        assertEquals(-M_PI + 0.1, r.project(-1e-15))
        assertEquals(0.0, full.project(0))
        assertEquals(M_PI, full.project(M_PI))
        assertEquals(M_PI, full.project(-M_PI))
    }

    fun testFromPointPair() {
        assertEquals(S1Interval.fromPointPair(-M_PI, M_PI), pi)
        assertEquals(S1Interval.fromPointPair(M_PI, -M_PI), pi)
        assertEquals(S1Interval.fromPointPair(mid34.hi, mid34.lo), mid34)
        assertEquals(S1Interval.fromPointPair(mid23.lo, mid23.hi), mid23)
    }

    fun testExpanded() {
        assertEquals(empty.expanded(1), empty)
        assertEquals(full.expanded(1), full)
        assertEquals(zero.expanded(1), S1Interval(-1.0, 1.0))
        assertEquals(mipi.expanded(0.01), S1Interval(M_PI - 0.01, -M_PI + 0.01))
        assertEquals(pi.expanded(27), full)
        assertEquals(pi.expanded(M_PI_2), quad23)
        assertEquals(pi2.expanded(M_PI_2), quad12)
        assertEquals(mipi2.expanded(M_PI_2), quad34)

        assertEquals(empty.expanded(-1), empty)
        assertEquals(full.expanded(-1), full)
        assertEquals(quad123.expanded(-27), empty)
        assertEquals(quad234.expanded(-27), empty)
        assertEquals(quad123.expanded(-M_PI_2), quad2)
        assertEquals(quad341.expanded(-M_PI_2), quad4)
        assertEquals(quad412.expanded(-M_PI_2), quad1)
    }

    fun testApproxEquals() {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        val kLo = 4 * DBL_EPSILON;  // < max_error default
        val kHi = 6 * DBL_EPSILON;  // > max_error default

        // Empty intervals.
        assertTrue(empty.approxEquals(empty))
        assertTrue(zero.approxEquals(empty) && empty.approxEquals(zero))
        assertTrue(pi.approxEquals(empty) && empty.approxEquals(pi))
        assertTrue(mipi.approxEquals(empty) && empty.approxEquals(mipi))
        assertFalse(empty.approxEquals(full))
        assertTrue(empty.approxEquals(S1Interval(1.0, 1 + 2 * kLo)))
        assertFalse(empty.approxEquals(S1Interval(1.0, 1 + 2 * kHi)))
        assertTrue(S1Interval(M_PI - kLo, -M_PI + kLo).approxEquals(empty))

        // Full intervals.
        assertTrue(full.approxEquals(full))
        assertFalse(full.approxEquals(empty))
        assertFalse(full.approxEquals(zero))
        assertFalse(full.approxEquals(pi))
        assertTrue(full.approxEquals(S1Interval(kLo, -kLo)))
        assertFalse(full.approxEquals(S1Interval(2 * kHi, 0.0)))
        assertTrue(S1Interval(-M_PI + kLo, M_PI - kLo).approxEquals(full))
        assertFalse(S1Interval(-M_PI, M_PI - 2 * kHi).approxEquals(full))

        // Singleton intervals.
        assertTrue(pi.approxEquals(pi) && mipi.approxEquals(pi))
        assertTrue(pi.approxEquals(S1Interval(M_PI - kLo, M_PI - kLo)))
        assertFalse(pi.approxEquals(S1Interval(M_PI - kHi, M_PI - kHi)))
        assertTrue(pi.approxEquals(S1Interval(M_PI - kLo, -M_PI + kLo)))
        assertFalse(pi.approxEquals(S1Interval(M_PI - kHi, -M_PI)))
        assertFalse(zero.approxEquals(pi))
        assertTrue(pi.union(mid12).union(zero).approxEquals(quad12))
        assertTrue(quad2.intersection(quad3).approxEquals(pi))
        assertTrue(quad3.intersection(quad2).approxEquals(pi))

        // Intervals whose corresponding endpoints are nearly the same but where the
        // endpoints are in opposite order (i.e., inverted intervals).
        assertFalse(S1Interval(0.0, kLo).approxEquals(S1Interval(kLo, 0.0)))
        assertFalse(S1Interval(M_PI - 0.5 * kLo, -M_PI + 0.5 * kLo).approxEquals(S1Interval(-M_PI + 0.5 * kLo, M_PI - 0.5 * kLo)))

        // Other intervals.
        assertTrue(S1Interval(1 - kLo, 2 + kLo).approxEquals(S1Interval(1, 2)))
        assertTrue(S1Interval(1 + kLo, 2 - kLo).approxEquals(S1Interval(1, 2)))
        assertTrue(S1Interval(2 - kLo, 1 + kLo).approxEquals(S1Interval(2, 1)))
        assertTrue(S1Interval(2 + kLo, 1 - kLo).approxEquals(S1Interval(2, 1)))
        assertFalse(S1Interval(1 - kHi, 2 + kLo).approxEquals(S1Interval(1, 2)))
        assertFalse(S1Interval(1 + kHi, 2 - kLo).approxEquals(S1Interval(1, 2)))
        assertFalse(S1Interval(2 - kHi, 1 + kLo).approxEquals(S1Interval(2, 1)))
        assertFalse(S1Interval(2 + kHi, 1 - kLo).approxEquals(S1Interval(2, 1)))
        assertFalse(S1Interval(1 - kLo, 2 + kHi).approxEquals(S1Interval(1, 2)))
        assertFalse(S1Interval(1 + kLo, 2 - kHi).approxEquals(S1Interval(1, 2)))
        assertFalse(S1Interval(2 - kLo, 1 + kHi).approxEquals(S1Interval(2, 1)))
        assertFalse(S1Interval(2 + kLo, 1 - kHi).approxEquals(S1Interval(2, 1)))
    }

    fun testGetDirectedHausdorffDistance() {
        assertEquals(0.0, empty.getDirectedHausdorffDistance(empty))
        assertEquals(0.0, empty.getDirectedHausdorffDistance(mid12))
        assertEquals(M_PI, mid12.getDirectedHausdorffDistance(empty))

        assertEquals(0.0, quad12.getDirectedHausdorffDistance(quad123))
        val interval = S1Interval(3.0, -3.0);  // an interval whose complement center is 0.
        assertEquals(3.0, S1Interval(-0.1, 0.2).getDirectedHausdorffDistance(interval))
        assertEquals(3.0 - 0.1, S1Interval(0.1, 0.2).getDirectedHausdorffDistance(interval))
        assertEquals(3.0 - 0.1, S1Interval(-0.2, -0.1).getDirectedHausdorffDistance(interval))
    }

}