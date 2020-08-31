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
package dilivia.s2

import com.google.common.geometry.GeometryTestCase
import com.google.common.geometry.S2.*
import kotlin.math.atan
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.tan

class S1ChordAngleTest : GeometryTestCase() {

    fun testDefaultConstructor() {
        // Check that the default constructor returns an angle of 0.
        val a = S1ChordAngle()
        assertEquals(S1ChordAngle.zero, a)
    }

    fun testTwoPointConstructor() {
        for (iter in 0 until 100) {
            val frame = randomFrame
            val x = frame[0]
            val y = frame[1]
            val z = frame[2]
            assertEquals(S1Angle.zero, S1ChordAngle.between(z, z).toAngle())
            assertEquals(M_PI, S1ChordAngle.between(S2Point.unaryMinus(z), z).radians(), 1e-7)
            assertEquals(M_PI_2, S1ChordAngle.between(x, z).radians(), 1e-15)
            val w = S2Point.normalize(S2Point.plus(y, z))
            assertEquals(M_PI_4, S1ChordAngle.between(w, z).radians(), 1e-15)
        }
    }

    fun testFromLength2() {
        assertEquals(0.0, S1ChordAngle.fromLength2(0).degrees(), 1e-15)
        assertEquals(60.0, S1ChordAngle.fromLength2(1).degrees(), 1e-14)
        assertEquals(90.0, S1ChordAngle.fromLength2(2).degrees(), 1e-13)
        assertEquals(180.0, S1ChordAngle.fromLength2(4).degrees(), 1e-15)
        assertEquals(180.0, S1ChordAngle.fromLength2(5).degrees(), 1e-15)
    }

    fun testZero() {
        assertEquals(S1Angle.zero, S1ChordAngle.zero.toAngle())
    }

    fun testRight() {
        assertEquals(90.0, S1ChordAngle.right.degrees(), 1e-13)
    }

    fun testStraight() {
        assertEquals(S1Angle.degrees(180.0), S1ChordAngle.straight.toAngle())
    }

    fun testInfinity() {
        assertTrue(S1ChordAngle.straight < S1ChordAngle.infinity)
        assertEquals(S1ChordAngle.infinity, S1ChordAngle.infinity)
        assertEquals(S1Angle.infinity, S1ChordAngle.infinity.toAngle())
    }

    fun testNegative() {
        assertTrue(S1ChordAngle.negative < S1ChordAngle.zero)
        assertEquals(S1ChordAngle.negative, S1ChordAngle.negative)
        assertTrue(S1ChordAngle.negative.toAngle() < S1Angle.zero)
    }

    fun testPredicates() {
        assertTrue(S1ChordAngle.zero.isZero())
        assertFalse(S1ChordAngle.zero.isNegative())
        assertFalse(S1ChordAngle.zero.isSpecial())
        assertFalse(S1ChordAngle.straight.isSpecial())
        assertTrue(S1ChordAngle.negative.isNegative())
        assertTrue(S1ChordAngle.negative.isSpecial())
        assertTrue(S1ChordAngle.infinity.isInfinity())
        assertTrue(S1ChordAngle.infinity.isSpecial())
    }

    fun testToFromS1Angle() {
        assertEquals(0.0, S1ChordAngle(S1Angle.zero).radians())
        assertEquals(4.0, S1ChordAngle(S1Angle.radians(M_PI)).length2)
        assertEquals(M_PI, S1ChordAngle(S1Angle.radians(M_PI)).radians())
        assertEquals(S1Angle.infinity, S1ChordAngle(S1Angle.infinity).toAngle())
        assertTrue(S1ChordAngle(S1Angle.radians(-1)).radians() < 0)
        assertEquals(1.0, S1ChordAngle(S1Angle.radians(1.0)).radians())
    }

    fun testSuccessor() {
        assertEquals(S1ChordAngle.zero, S1ChordAngle.negative.successor())
        assertEquals(S1ChordAngle.infinity, S1ChordAngle.straight.successor())
        assertEquals(S1ChordAngle.infinity, S1ChordAngle.infinity.successor())
        var x = S1ChordAngle.negative
        for (i in 0 until 10) {
            assertTrue(x < x.successor())
            x = x.successor()
        }
    }

    fun testPredecessor() {
        assertEquals(S1ChordAngle.straight, S1ChordAngle.infinity.predecessor())
        assertEquals(S1ChordAngle.negative, S1ChordAngle.zero.predecessor())
        assertEquals(S1ChordAngle.negative, S1ChordAngle.negative.predecessor())
        var x = S1ChordAngle.infinity
        for (i in 0 until 10) {
            assertTrue(x > x.predecessor())
            x = x.predecessor()
        }
    }

    fun testArithmetic() {
        val zero = S1ChordAngle.zero
        val degree30 = S1ChordAngle.degrees(30)
        val degree60 = S1ChordAngle.degrees(60)
        val degree90 = S1ChordAngle.degrees(90)
        val degree120 = S1ChordAngle.degrees(120)
        val degree180 = S1ChordAngle.straight
        assertEquals(0.0, (zero + zero).degrees(), 1e-15)
        assertEquals(0.0, (zero - zero).degrees(), 1e-15)
        assertEquals(0.0, (degree60 - degree60).degrees(), 1e-15)
        assertEquals(0.0, (degree180 - degree180).degrees(), 1e-15)
        assertEquals(0.0, (zero - degree60).degrees(), 1e-15)
        assertEquals(0.0, (degree30 - degree90).degrees(), 1e-15)
        assertEquals(60.0, (degree60 + zero).degrees(), 1e-14)
        assertEquals(60.0, (degree60 - zero).degrees(), 1e-14)
        assertEquals(60.0, (zero + degree60).degrees(), 1e-14)
        assertEquals(90.0, (degree30 + degree60).degrees(), 1e-13)
        assertEquals(90.0, (degree60 + degree30).degrees(), 1e-13)
        assertEquals(60.0, (degree90 - degree30).degrees(), 1e-14)
        assertEquals(30.0, (degree90 - degree60).degrees(), 1e-14)
        assertEquals(180.0, (degree180 + zero).degrees(), 1e-15)
        assertEquals(180.0, (degree180 - zero).degrees(), 1e-15)
        assertEquals(180.0, (degree90 + degree90).degrees(), 1e-15)
        assertEquals(180.0, (degree120 + degree90).degrees(), 1e-15)
        assertEquals(180.0, (degree120 + degree120).degrees(), 1e-15)
        assertEquals(180.0, (degree30 + degree180).degrees(), 1e-15)
        assertEquals(180.0, (degree180 + degree180).degrees(), 1e-15)
    }

    fun testTrigonometry() {
        val kIters = 20
        for (iter in 0..kIters) {
            val radians = M_PI * iter / kIters
            val angle = S1ChordAngle(S1Angle.radians(radians))
            assertEquals(sin(radians), sin(angle), 1e-15)
            assertEquals(cos(radians), cos(angle), 1e-15)
            // Since the tan(x) is unbounded near Pi/4, we map the result back to an
            // angle before comparing.  (The assertion is that the result is equal to
            // the tangent of a nearby angle.)
            assertEquals(atan(tan(radians)), atan(tan(angle)), 1e-15)
        }

        // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
        val angle90 = S1ChordAngle.fromLength2(2)
        val angle180 = S1ChordAngle.fromLength2(4)
        assertEquals(1.0, sin(angle90))
        assertEquals(0.0, cos(angle90))
        assertEquals(Double.POSITIVE_INFINITY, tan(angle90))
        assertEquals(0.0, sin(angle180))
        assertEquals(-1.0, cos(angle180))
        assertEquals(0.0, tan(angle180))
    }

    fun testPlusError() {
        assertEquals(S1ChordAngle.negative, S1ChordAngle.negative.plusError(5))
        assertEquals(S1ChordAngle.infinity, S1ChordAngle.infinity.plusError(-5))
        assertEquals(S1ChordAngle.straight, S1ChordAngle.straight.plusError(5))
        assertEquals(S1ChordAngle.zero, S1ChordAngle.zero.plusError(-5))
        assertEquals(S1ChordAngle.fromLength2(1.25), S1ChordAngle.fromLength2(1).plusError(0.25))
        assertEquals(S1ChordAngle.fromLength2(0.75), S1ChordAngle.fromLength2(1).plusError(-0.25))
    }
/*
    fun testGetS2PointConstructorMaxError() {
        // Check that the error bound returned by GetS2PointConstructorMaxError() is
        // large enough.
        for (iter in 0 until 100000) {
            rand.setSeed(iter.toLong()) // Easier to reproduce a specific case.
            val x = randomPoint()
            var y = randomPoint()
            if (Random.Default.nextInt(0, 10) == 0) {
                // Occasionally test a point pair that is nearly identical or antipodal.
                val r = S1Angle . radians (1e-15 * rand.nextDouble())
                y = S2EdgeDistances.interpolateAtDistance(r, x, y)
                if (rand.nextBoolean()) y = -y
            }
            val dist = S1ChordAngle (x, y)
            val error = dist . getS2PointConstructorMaxError ()
            assertTrue(S2Predicates.compareDistance(x, y, dist.plusError(error)) <= 0.0)
            assertTrue(S2Predicates.compareDistance(x, y, dist.plusError(-error)) >= 0.0)
        }
    }
    */


}