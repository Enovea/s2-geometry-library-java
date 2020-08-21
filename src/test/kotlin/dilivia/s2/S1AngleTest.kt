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

import com.google.common.geometry.S2.*
import com.google.common.geometry.S2LatLng
import dilivia.s2.S1Angle.Companion.degrees
import dilivia.s2.S1Angle.Companion.e5
import dilivia.s2.S1Angle.Companion.e6
import dilivia.s2.S1Angle.Companion.e7
import dilivia.s2.S1Angle.Companion.radians
import junit.framework.TestCase
import kotlin.random.Random

@Strictfp
class S1AngleTest : TestCase() {

    fun testBasic() {
        // Check that the conversion between Pi radians and 180 degrees is exact.
        assertEquals(radians(Math.PI).radians, Math.PI)
        assertEquals(radians(Math.PI).degrees(), 180.0)
        assertEquals(degrees(180.0).radians, Math.PI)
        assertEquals(degrees(180.0).degrees(), 180.0)
        assertEquals(radians(Math.PI / 2).degrees(), 90.0)

        // Check negative angles.
        assertEquals(radians(-Math.PI / 2).degrees(), -90.0)
        assertEquals(degrees(-45.0).radians, -Math.PI / 4)

        // Check that E5/E6/E7 representations work as expected.
        assertEquals(e5(2000000), degrees(20.0))
        assertEquals(e6(-60000000), degrees(-60.0))
        assertEquals(e7(750000000), degrees(75.0))
        assertEquals(degrees(12.34567).e5(), 1234567)
        assertEquals(degrees(12.345678).e6(), 12345678)
        assertEquals(degrees(-12.3456789).e7(), -123456789)
    }

    fun testDefaultConstructor() {
        // Check that the default constructor returns an angle of 0.
        assertEquals(0.0, S1Angle().radians)
    }

    fun testInfinity() {
        assertTrue(radians(1e30) < S1Angle.infinity)
        assertTrue(-S1Angle.infinity < S1Angle.zero)
        assertEquals(S1Angle.infinity, S1Angle.infinity)
    }

    fun testZero() {
        assertEquals(radians(0), S1Angle.zero)
    }

    fun testPiRadiansExactly180Degrees() {
        // Check that the conversion between Pi radians and 180 degrees is exact.
        assertEquals(M_PI, radians(M_PI).radians)
        assertEquals(180.0, radians(M_PI).degrees())
        assertEquals(M_PI, degrees(180).radians)
        assertEquals(180.0, degrees(180).degrees())

        assertEquals(90.0, radians(M_PI_2).degrees())

        // Check negative angles.
        assertEquals(-90.0, radians(-M_PI_2).degrees())
        assertEquals(-M_PI_4, degrees(-45).radians)
    }

    fun testE5E6E7Representations() {
        // Check that E5/E6/E7 representations work as expected.
        assertEquals(degrees(-45).radians, e5(-4500000).radians, DBL_EPSILON)
        assertEquals(degrees(-60).radians, e6(-60000000).radians, DBL_EPSILON)
        assertEquals(degrees(75).radians, e7(750000000).radians, DBL_EPSILON)
        assertEquals(-17256123, degrees(-172.56123).e5())
        assertEquals(12345678, degrees(12.345678).e6())
        assertEquals(-123456789, degrees(-12.3456789).e7())
    }

    fun testNormalizeCorrectlyCanonicalizesAngles() {
        assertEquals(0.0, degrees(360.0).normalized().degrees())
        assertEquals(-90.0, degrees(-90.0).normalized().degrees())
        assertEquals(180.0, degrees(-180.0).normalized().degrees())
        assertEquals(180.0, degrees(180.0).normalized().degrees())
        assertEquals(180.0, degrees(540.0).normalized().degrees())
        assertEquals(90.0, degrees(-270.0).normalized().degrees())
    }

    fun testArithmeticOperationsOnAngles() {
        assertEquals(0.3, radians(-0.3).abs().radians, DBL_EPSILON)
        assertEquals(-0.1, (-radians(0.1)).radians, DBL_EPSILON)
        assertEquals(0.4, (radians(0.1) + radians(0.3)).radians, DBL_EPSILON)
        assertEquals(-0.2, (radians(0.1) - radians(0.3)).radians, DBL_EPSILON)
        assertEquals(0.6, (2.0 * radians(0.3)).radians, DBL_EPSILON)
        assertEquals(0.6, (radians(0.3) * 2.0).radians, DBL_EPSILON)
        assertEquals(0.15, (radians(0.3) / 2.0).radians, DBL_EPSILON)
        assertEquals(0.5, (radians(0.3) / radians(0.6)).radians, DBL_EPSILON)

        var tmp = radians(1.0)
        tmp += radians(0.5)
        assertEquals(1.5, tmp.radians)
        tmp -= radians(1.0)
        assertEquals(0.5, tmp.radians)
        tmp *= 5.0
        assertEquals(2.5, tmp.radians)
        tmp /= 2.0
        assertEquals(1.25, tmp.radians)
    }

    fun testTrigonometry() {
        // Spot check a few angles to ensure that the correct function is called.
        assertEquals(1.0, cos(degrees(0)), DBL_EPSILON)
        assertEquals(1.0, sin(degrees(90)), DBL_EPSILON)
        assertEquals(1.0, tan(degrees(45)), DBL_EPSILON)
    }

    fun testConstructorsThatMeasureAngles() {
        assertEquals(M_PI_2, S1Angle(S2Point(1, 0, 0), S2Point(0, 0, 2)).radians)
        assertEquals(0.0, S1Angle(S2Point(1, 0, 0), S2Point(1, 0, 0)).radians)
        assertEquals(50.0, S1Angle(S2LatLng.fromDegrees(20, 20), S2LatLng.fromDegrees(70, 20)).degrees(), 1e-13)
    }

    fun testTextFormatting() {
        assertEquals("180.0000000d", degrees(180.0).toString())
    }

    fun testDegreesVsE6() {
// The current implementation guarantees exact conversions between
// Degrees() and E6() when the Degrees() argument is an integer.
        for (i in 0..180) {
            assertEquals(degrees(i), e6(1000000L * i))
        }
    }

    fun testDegreesVsE7() {
// The current implementation guarantees exact conversions between
// Degrees() and E7() when the Degrees() argument is an integer.
        for (i in 0..180) {
            assertEquals(degrees(i), e7(10000000L * i))
        }
    }

    fun testE6VsE7() {
// The current implementation guarantees exact conversions between
// E6() and E7() when the E6() argument is an integer.
        for (iter in 0..1000) {
            val i = Random.Default.nextLong(0, 180000000)
            assertEquals(e6(i), e7(10 * i))
        }
    }

    fun testDegreesVsRadians() {
// The current implementation guarantees certain exact conversions between
// degrees and radians (see the header file for details).
        for (k in -8..8) {
        assertEquals(degrees(45 * k), radians(k * M_PI / 4))
            assertEquals(45.0 * k, degrees(45 * k).degrees())
        }
        for (k in 0..30) {
        val n = 1 shl k
        assertEquals(degrees(180.0/ n), radians(M_PI / n))
            assertEquals(degrees(60.0/ n), radians(M_PI / (3.0* n)))
            assertEquals(degrees(36.0/ n), radians(M_PI / (5.0* n)))
            assertEquals(degrees(20.0/ n), radians(M_PI / (9.0* n)))
            assertEquals(degrees(4.0/ n), radians(M_PI / (45.0* n)))
        }
        // We also spot check a couple of non-identities.
        assertTrue(degrees(3) != radians(M_PI / 60))
        assertTrue(60.0 != degrees(60).degrees())
    }

}