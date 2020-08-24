/*
 * Copyright 2005 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
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
class S2LatLngTest : GeometryTestCase() {

    fun testBasic() {
        val llRad = S2LatLng.fromRadians(M_PI_4, M_PI_2)
        assertEquals(M_PI_4, llRad.lat().radians)
        assertEquals(M_PI_2, llRad.lng().radians)
        assertTrue(llRad.isValid)
        val llDeg = S2LatLng.fromDegrees(45, 90)
        assertEquals(llRad, llDeg)
        assertTrue(llDeg.isValid)
        assertFalse(S2LatLng.fromDegrees(-91, 0).isValid)
        assertFalse(S2LatLng.fromDegrees(0, 181).isValid)

        var bad = S2LatLng.fromDegrees(120, 200)
        assertFalse(bad.isValid)
        var better = bad.normalized()
        assertTrue(better.isValid)
        assertEquals(S1Angle.degrees(90), better.lat())
        assertEquals(S1Angle.degrees(-160).radians, better.lng().radians)

        bad = S2LatLng.fromDegrees(-100, -360)
        assertFalse(bad.isValid)
        better = bad.normalized()
        assertTrue(better.isValid)
        assertEquals(S1Angle.degrees(-90), better.lat())
        assertEquals(0.0, better.lng().radians, 1e-15)

        assertTrue((S2LatLng.fromDegrees(10, 20) + S2LatLng.fromDegrees(20, 30)).approxEquals(S2LatLng.fromDegrees(30, 50)))
        assertTrue((S2LatLng.fromDegrees(10, 20) - S2LatLng.fromDegrees(20, 30)).approxEquals(S2LatLng.fromDegrees(-10, -10)))
        assertTrue((0.5 * S2LatLng.fromDegrees(10, 20)).approxEquals(S2LatLng.fromDegrees(5, 10)))

        // Check that Invalid() returns an invalid point.
        assertFalse(S2LatLng.invalid.isValid)

        // Check that the default constructor sets latitude and longitude to 0.
        val center = S2LatLng.center
        assertTrue(center.isValid)
        assertEquals(0.0, center.lat().radians)
        assertEquals(0.0, center.lng().radians)
    }

    fun testConversion() {
        // Test special cases: poles, "date line"
        assertEquals(90.0, S2LatLng.fromPoint(S2LatLng.fromDegrees(90.0, 65.0).toPoint()).lat().degrees())
        assertEquals(-M_PI_2, S2LatLng.fromPoint(S2LatLng.fromRadians(-M_PI_2, 1.0).toPoint()).lat().radians)
        assertEquals(180.0, abs(S2LatLng.fromPoint(S2LatLng.fromDegrees(12.2, 180.0).toPoint()).lng().degrees()))
        assertEquals(M_PI, abs(S2LatLng.fromPoint(S2LatLng.fromRadians(0.1, -M_PI).toPoint()).lng().radians))

        // Test a bunch of random points.
        for (i in 0 until 100000) {
            val p = randomPoint()
            assertTrue(p.approxEquals(S2LatLng.fromPoint(p).toPoint()))
        }
    }

    fun testDistance() {
        assertEquals(0.0, S2LatLng.fromDegrees(90, 0).getDistance(S2LatLng.fromDegrees(90, 0)).radians)
        assertEquals(77.0, S2LatLng.fromDegrees(-37, 25).getDistance(S2LatLng.fromDegrees(-66, -155)).degrees(), 1e-13)
        assertEquals(115.0, S2LatLng.fromDegrees(0, 165).getDistance(S2LatLng.fromDegrees(0, -80)).degrees(), 1e-13)
        assertEquals(180.0, S2LatLng.fromDegrees(47, -127).getDistance(S2LatLng.fromDegrees(-47, 53)).degrees(), 2e-6)
    }

}