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

import dilivia.s2.S2.M_PI
import dilivia.s2.coords.S2Coords
import dilivia.s2.region.S2Cell
import mu.KotlinLogging
import java.lang.Math.pow
import kotlin.math.IEEErem
import kotlin.math.acos
import kotlin.math.asin
import kotlin.math.pow

class S2PointTest : S2GeometryTestCase() {

    fun testFrames() {
        val z = S2Point(0.2, 0.5, -3.3).normalize()
        val m = S2Point.getFrame(z)
        assertTrue(S2Point.approxEquals(m.col(2).toS2Point(), z))
        assertTrue(S2Point.isUnitLength(m.col(0).toS2Point()))
        assertTrue(S2Point.isUnitLength(m.col(1).toS2Point()))
        assertEquals(m.det(), 1.0, 1e-15)

        assertTrue(S2Point.approxEquals(S2Point.toFrame(m, m.col(0).toS2Point()), S2Point(1, 0, 0)))
        assertTrue(S2Point.approxEquals(S2Point.toFrame(m, m.col(1).toS2Point()), S2Point(0, 1, 0)))
        assertTrue(S2Point.approxEquals(S2Point.toFrame(m, m.col(2).toS2Point()), S2Point(0, 0, 1)))

        assertTrue(S2Point.approxEquals(S2Point.fromFrame(m, S2Point(1, 0, 0)), m.col(0).toS2Point()))
        assertTrue(S2Point.approxEquals(S2Point.fromFrame(m, S2Point(0, 1, 0)), m.col(1).toS2Point()))
        assertTrue(S2Point.approxEquals(S2Point.fromFrame(m, S2Point(0, 0, 1)), m.col(2).toS2Point()))
    }

    fun testRotate() {
        repeat(1000) {
            val axis = S2Random.randomPoint()
            val target = S2Random.randomPoint()
            // Choose a distance whose logarithm is uniformly distributed.
            var distance = M_PI * 1e-15.pow(S2Random.randomDouble())
            // Sometimes choose points near the far side of the axis.
            if (S2Random.oneIn(5)) distance = M_PI - distance
            val p = S2EdgeDistances.interpolateAtDistance(S1Angle.radians(distance), axis, target)
            // Choose the rotation angle.
            var angle = 2 * M_PI * pow(1e-15, S2Random.randomDouble())
            if (S2Random.oneIn(3)) angle = -angle
            if (S2Random.oneIn(10)) angle = 0.0
            testRotate(p, axis, S1Angle.radians(angle))
        }
    }


    fun testOriginTest() {
        // To minimize the number of expensive Sign() calculations,
        // S2::Origin() should not be nearly collinear with any commonly used edges.
        // Two important categories of such edges are:
        //
        //  - edges along a line of longitude (reasonably common geographically)
        //  - S2Cell edges (used extensively when computing S2Cell coverings)
        //
        // This implies that the origin:
        //
        //  - should not be too close to either pole (since all lines of longitude
        //    converge at the poles)
        //  - should not be colinear with edges of any S2Cell except for very small
        //    ones (which are used less frequently)
        //
        // The point chosen below is about 66km from the north pole towards the East
        // Siberian Sea.  The purpose of the STtoUV(2/3) calculation is to keep the
        // origin as far away as possible from the longitudinal edges of large
        // S2Cells.  (The line of longitude through the chosen point is always 1/3
        // or 2/3 of the way across any S2Cell with longitudinal edges that it
        // passes through.)

      assertEquals(S2Point.origin(), S2Point(-0.01, 0.01 * S2Coords.stToUV(2.0/3.0), 1.0).normalize())

        // Check that the origin is not too close to either pole.  (We don't use
        // S2Earth because we don't want to depend on that package.)
        val distance_km = acos (S2Point.origin().z()) * kEarthRadiusKm
        assertTrue(distance_km >= 50.0)
        logger.info { "\nS2::Origin() coordinates: ${S2LatLng.fromPoint(S2Point.origin())}, distance from pole: $distance_km km" }

        // Check that S2::Origin() is not collinear with the edges of any large
        // S2Cell.  We do this is two parts.  For S2Cells that belong to either
        // polar face, we simply need to check that S2::Origin() is not nearly
        // collinear with any edge of any cell that contains it (except for small
        // cells < 3 meters across).
        assertTrue(getMinExpensiveLevel(S2Point.origin()) >= 22)

        // For S2Cells that belong to the four non-polar faces, only longitudinal
        // edges can possibly be colinear with S2::Origin().  We check these edges
        // by projecting S2::Origin() onto the equator, and then testing all S2Cells
        // that contain this point to make sure that none of their edges are nearly
        // colinear with S2::Origin() (except for small cells < 3 meters across).
        val equator_point = S2Point(S2Point.origin().x(), S2Point.origin().y(), 0.0)
      assertTrue(getMinExpensiveLevel(equator_point) >= 22)
    }

    companion object {

      val logger = KotlinLogging.logger {  }

        fun testRotate(p: S2Point, axis: S2Point, angle: S1Angle) {
            val result = S2Point.rotate(p, axis, angle)

            // "result" should be unit length.
            assertTrue(S2Point.isUnitLength(result))

            // "result" and "p" should be the same distance from "axis".
            val kMaxPositionError = 1e-15
            assertTrue((S1Angle(result, axis) - S1Angle(p, axis)).abs().radians <= kMaxPositionError)

            // Check that the rotation angle is correct.  We allow a fixed error in the
            // *position* of the result, so we need to convert this into a rotation
            // angle.  The allowable error can be very large as "p" approaches "axis".
            val axis_distance = p.crossProd(axis).norm()
            val max_rotation_error =
                    if (axis_distance < kMaxPositionError) {
                        2 * M_PI
                    } else {
                        asin(kMaxPositionError / axis_distance)
                    }
            val actual_rotation = S2Measures.turnAngle(p, axis, result) + M_PI
            val rotation_error = (angle.radians - actual_rotation).IEEErem(2 * M_PI)
            assertTrue(rotation_error <= max_rotation_error)
        }

        // Given a point P, return the minimum level at which an edge of some S2Cell
        // parent of P is nearly collinear with S2::Origin().  This is the minimum
        // level for which Sign() may need to resort to expensive calculations in
        // order to determine which side of an edge the origin lies on.
        fun getMinExpensiveLevel(p: S2Point): Int {
            val id = S2CellId.fromPoint(p)
            for (level in 0..S2CellId.kMaxLevel) {
                val cell = S2Cell(id.parent(level))
                for (k in 0..3) {
                    val a = cell.getVertex(k)
                    val b = cell.getVertex(k + 1)
                    if (S2Predicates.triageSign(a, b, S2Point.origin(), a.crossProd(b)) == 0) {
                        return level
                    }
                }
            }
            return S2CellId.kMaxLevel + 1
        }

    }
}
