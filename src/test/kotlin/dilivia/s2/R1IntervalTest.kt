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
import com.google.common.geometry.S2.DBL_EPSILON
import dilivia.s2.R1Interval.Companion.empty
import dilivia.s2.R1Interval.Companion.fromPointPair


@Strictfp
class R1IntervalTest : GeometryTestCase() {
    /**
     * Test all of the interval operations on the given pair of intervals.
     * "expected_relation" is a sequence of "T" and "F" characters corresponding
     * to the expected results of contains(), interiorContains(), Intersects(),
     * and InteriorIntersects() respectively.
     */
    private fun testIntervalOps(x: R1Interval, y: R1Interval?, expectedRelation: String) {
        assertEquals(x.contains(y!!), expectedRelation[0] == 'T')
        assertEquals(x.interiorContains(y), expectedRelation[1] == 'T')
        assertEquals(x.intersects(y), expectedRelation[2] == 'T')
        assertEquals(x.interiorIntersects(y), expectedRelation[3] == 'T')
        assertEquals(x.contains(y), x.union(y) == x)
        assertEquals(x.intersects(y), !x.intersection(y).isEmpty)

        val z = x.addInterval(y)
        assertEquals(x.union(y), z)
    }

    fun testBasic() {
        // Constructors and accessors.
        val unit = R1Interval(0, 1)
        val negunit = R1Interval(-1, 0)
        assertEquals(unit.lo, 0.0)
        assertEquals(unit.hi, 1.0)
        assertEquals(negunit.lo, -1.0)
        assertEquals(negunit.hi, 0.0)

        // Mutable
        val ten = MutableR1Interval(0, 0)
        ten.hi = 10.0
        assertEquals(R1Interval(0, 10), ten);
        ten[0] = -10.0
        assertEquals(R1Interval(-10, 10), ten)
        ten[1] = 0.0
        assertEquals(R2Vector(-10.0, 0.0), ten.bounds)
        ten.bounds = R2Vector(0, 10)
        assertEquals(R1Interval(0, 10), ten);

        // is_empty()
        val half = R1Interval(0.5, 0.5)
        assertTrue(unit.isNotEmpty)
        assertFalse(unit.isEmpty)
        assertTrue(half.isNotEmpty)
        assertFalse(half.isEmpty)
        val empty = empty()
        assertTrue(empty.isEmpty)
        assertFalse(empty.isNotEmpty)

        // GetCenter(), GetLength()
        assertEquals(unit.center, 0.5)
        assertEquals(half.center, 0.5)
        assertEquals(negunit.length, 1.0)
        assertEquals(half.length, 0.0)
        assertTrue(empty.length < 0)

        // == and !=
        assertTrue(empty == empty)
        assertTrue(unit == unit)
        assertTrue(unit != empty)
        assertTrue(R1Interval(1, 2) != R1Interval(1, 3))

        // Check that the default R1Interval is identical to Empty().
        val defaultEmpty = R1Interval()
        assertTrue(defaultEmpty.isEmpty)
        assertEquals(empty.lo, defaultEmpty.lo)
        assertEquals(empty.hi, defaultEmpty.hi)

        // GetCenter(), GetLength()
        assertEquals(unit.center, 0.5);
        assertEquals(half.center, 0.5);
        assertEquals(negunit.length, 1.0);
        assertEquals(half.length, 0.0);
        assertTrue(empty.length < 0);

        // contains(double), interiorContains(double)
        assertTrue(unit.contains(0.5))
        assertTrue(unit.interiorContains(0.5))
        assertTrue(unit.contains(0.0))
        assertTrue(!unit.interiorContains(0.0))
        assertTrue(unit.contains(1.0))
        assertTrue(!unit.interiorContains(1.0))

        // contains(R1Interval), interiorContains(R1Interval)
        // Intersects(R1Interval), InteriorIntersects(R1Interval)
        testIntervalOps(empty, empty, "TTFF")
        testIntervalOps(empty, unit, "FFFF")
        testIntervalOps(unit, half, "TTTT")
        testIntervalOps(unit, unit, "TFTT")
        testIntervalOps(unit, empty, "TTFF")
        testIntervalOps(unit, negunit, "FFTF")
        testIntervalOps(unit, R1Interval(0.0, 0.5), "TFTT")
        testIntervalOps(half, R1Interval(0.0, 0.5), "FFTF")

        // addPoint()
        var r: R1Interval = empty.addPoint(5.0)
        assertTrue(r.lo == 5.0 && r.hi == 5.0)
        r = r.addPoint(-1.0)
        assertTrue(r.lo == -1.0 && r.hi == 5.0)
        r = r.addPoint(0.0)
        assertTrue(r.lo == -1.0 && r.hi == 5.0)

        // Project()
        assertEquals(0.3, R1Interval(0.1, 0.4).project(0.3));
        assertEquals(0.1, R1Interval(0.1, 0.4).project(-7.0));
        assertEquals(0.4, R1Interval(0.1, 0.4).project(0.6));

        // fromPointPair()
        assertEquals(fromPointPair(4.0, 4.0), R1Interval(4, 4))
        assertEquals(fromPointPair(-1.0, -2.0), R1Interval(-2, -1))
        assertEquals(fromPointPair(-5.0, 3.0), R1Interval(-5, 3))

        // expanded()
        assertEquals(empty.expanded(0.45), empty)
        assertEquals(unit.expanded(0.5), R1Interval(-0.5, 1.5))
        assertEquals(R1Interval(0.5, 0.5), unit.expanded(-0.5));
        assertTrue(unit.expanded(-0.51).isEmpty);
        assertTrue(unit.expanded(-0.51).expanded(0.51).isEmpty);

        // union(), intersection()
        assertTrue(R1Interval(99, 100).union(empty) == R1Interval(99, 100))
        assertTrue(empty.union(R1Interval(99, 100)) == R1Interval(99, 100))
        assertTrue(R1Interval(5, 3).union(R1Interval(0, -2)).isEmpty)
        assertTrue(R1Interval(0, -2).union(R1Interval(5, 3)).isEmpty)
        assertTrue(unit.union(unit) == unit)
        assertTrue(unit.union(negunit) == R1Interval(-1, 1))
        assertTrue(negunit.union(unit) == R1Interval(-1, 1))
        assertTrue(half.union(unit) == unit)
        assertTrue(unit.intersection(half) == half)
        assertTrue(unit.intersection(negunit) == R1Interval(0, 0))
        assertTrue(negunit.intersection(half).isEmpty)
        assertTrue(unit.intersection(empty).isEmpty)
        assertTrue(empty.intersection(unit).isEmpty)
    }

    fun testApproxEquals() {
        // Choose two values kLo and kHi such that it's okay to shift an endpoint by
        // kLo (i.e., the resulting interval is equivalent) but not by kHi.
        val kLo = 4 * DBL_EPSILON  // < max_error default
        val kHi = 6 * DBL_EPSILON  // > max_error default

        // Empty intervals.
        val empty = R1Interval.empty()
        assertTrue(empty.approxEquals(empty))
        assertTrue(R1Interval(0, 0).approxEquals(empty))
        assertTrue(empty.approxEquals(R1Interval(0, 0)))
        assertTrue(R1Interval(1, 1).approxEquals(empty))
        assertTrue(empty.approxEquals(R1Interval(1, 1)))
        assertFalse(empty.approxEquals(R1Interval(0, 1)))
        assertTrue(empty.approxEquals(R1Interval(1.0, 1 + 2*kLo)))
        assertFalse(empty.approxEquals(R1Interval(1.0, 1 + 2*kHi)))

        // Singleton intervals.
        assertTrue(R1Interval(1, 1).approxEquals(R1Interval(1, 1)))
        assertTrue(R1Interval(1, 1).approxEquals(R1Interval(1 - kLo, 1 - kLo)))
        assertTrue(R1Interval(1, 1).approxEquals(R1Interval(1 + kLo, 1 + kLo)))
        assertFalse(R1Interval(1, 1).approxEquals(R1Interval(1.0 - kHi, 1.0)))
        assertFalse(R1Interval(1, 1).approxEquals(R1Interval(1.0, 1 + kHi)))
        assertTrue(R1Interval(1, 1).approxEquals(R1Interval(1 - kLo, 1 + kLo)))
        assertFalse(R1Interval(0, 0).approxEquals(R1Interval(1, 1)))

        // Other intervals.
        assertTrue(R1Interval(1 - kLo, 2 + kLo).approxEquals(R1Interval(1, 2)))
        assertTrue(R1Interval(1 + kLo, 2 - kLo).approxEquals(R1Interval(1, 2)))
        assertFalse(R1Interval(1 - kHi, 2 + kLo).approxEquals(R1Interval(1, 2)))
        assertFalse(R1Interval(1 + kHi, 2 - kLo).approxEquals(R1Interval(1, 2)))
        assertFalse(R1Interval(1 - kLo, 2 + kHi).approxEquals(R1Interval(1, 2)))
        assertFalse(R1Interval(1 + kLo, 2 - kHi).approxEquals(R1Interval(1, 2)))
    }

    
}