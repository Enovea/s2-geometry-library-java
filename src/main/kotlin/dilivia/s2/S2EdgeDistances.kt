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

import com.google.common.geometry.S2
import kotlin.math.*

object S2EdgeDistances {

    /////////////////////////////////////////////////////////////////////////////
    ///////////////         (point along edge) functions          ///////////////

    // Given a point X and an edge AB, returns the distance ratio AX / (AX + BX).
    // If X happens to be on the line segment AB, this is the fraction "t" such
    // that X == Interpolate(t, A, B).  Requires that A and B are distinct.
    fun getDistanceFraction(x: S2Point, a: S2Point, b: S2Point): Double {
        assert(a <= b)
        val d0 = x.angle(a)
        val d1 = x.angle(b)
        return d0 / (d0 + d1)
    }

    // Returns the point X along the line segment AB whose distance from A is the
    // given fraction "t" of the distance AB.  Does NOT require that "t" be
    // between 0 and 1.  Note that all distances are measured on the surface of
    // the sphere, so this is more complicated than just computing (1-t)*a + t*b
    // and normalizing the result.
    fun interpolate(t: Double, a: S2Point, b: S2Point): S2Point {
        if (t == 0.0) return a
        if (t == 1.0) return b
        val ab = S1Angle(a, b)
        return interpolateAtDistance(t * ab, a, b)
    }

    // Like Interpolate(), except that the parameter "ax" represents the desired
    // distance from A to the result X rather than a fraction between 0 and 1.
    fun interpolateAtDistance(ax_angle: S1Angle, a: S2Point, b: S2Point): S2Point {
        val ax = ax_angle.radians

        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)

        // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
        // result is always perpendicular to A, even if A=B or A=-B, but it is not
        // necessarily unit length.  (We effectively normalize it below.)
        val normal = S2.robustCrossProd(a, b)
        val tangent = normal.crossProd(a)
        assert(tangent != S2Point(0, 0, 0))

        // Now compute the appropriate linear combination of A and "tangent".  With
        // infinite precision the result would always be unit length, but we
        // normalize it anyway to ensure that the error is within acceptable bounds.
        // (Otherwise errors can build up when the result of one interpolation is
        // fed into another interpolation.)
        return (cos(ax) * a + (sin(ax) / tangent.norm()) * tangent).normalize()
    }

    /////////////////////////////////////////////////////////////////////////////
    ///////////////            (point, edge) functions            ///////////////

    // Returns the minimum distance from X to any point on the edge AB.  All
    // arguments should be unit length.  The result is very accurate for small
    // distances but may have some numerical error if the distance is large
    // (approximately Pi/2 or greater).  The case A == B is handled correctly.
    //
    // If you want to compare a distance against a fixed threshold, e.g.
    //    if (S2::GetDistance(x, a, b) < limit)
    // then it is significantly faster to use UpdateMinDistance() below.
    fun getDistance(x: S2Point, a: S2Point, b: S2Point): S1Angle {
        val minDist = MutableS1ChordAngle(Double.MAX_VALUE)
        alwaysUpdateMinDistance(x, a, b, minDist, true)
        return minDist.toAngle()
    }

    // Returns true if the distance from X to the edge AB is less than "limit".
    // (Specify limit.Successor() for "less than or equal to".)  This method is
    // significantly faster than GetDistance().  If you want to compare against a
    // fixed S1Angle, you should convert it to an S1ChordAngle once and save the
    // value, since this step is relatively expensive.
    //
    // See s2pred::CompareEdgeDistance() for an exact version of this predicate.
    fun isDistanceLess(x: S2Point, a: S2Point, b: S2Point, limit: S1ChordAngle): Boolean = TODO()

    // If the distance from X to the edge AB is less than "min_dist", this
    // method updates "min_dist" and returns true.  Otherwise it returns false.
    // The case A == B is handled correctly.
    //
    // Use this method when you want to compute many distances and keep track of
    // the minimum.  It is significantly faster than using GetDistance(),
    // because (1) using S1ChordAngle is much faster than S1Angle, and (2) it
    // can save a lot of work by not actually computing the distance when it is
    // obviously larger than the current minimum.
    fun updateMinDistance(x: S2Point, a: S2Point, b: S2Point, minDistance: MutableS1ChordAngle): Boolean = alwaysUpdateMinDistance(x, a, b, minDistance, false)

    // If the maximum distance from X to the edge AB is greater than "max_dist",
    // this method updates "max_dist" and returns true.  Otherwise it returns false.
    // The case A == B is handled correctly.
    fun updateMaxDistance(x: S2Point, a: S2Point, b: S2Point, maxDistance: MutableS1ChordAngle): Boolean {
        val dist = maxOf(S1ChordAngle.between(x, a), S1ChordAngle.between(x, b)).toMutable()
        if (dist > S1ChordAngle.right) {
            alwaysUpdateMinDistance(-x, a, b, dist, true)
            dist.length2 = (S1ChordAngle.straight - dist).length2
        }
        if (maxDistance < dist) {
            maxDistance.length2 = dist.length2;
            return true;
        }

        return false;
    }

    // This function computes the distance from a point X to a line segment AB.
    // If the distance is less than "min_dist" or "always_update" is true, it
    // updates "min_dist" and returns true.  Otherwise it returns false.
    //
    // The "Always" in the function name refers to the template argument, i.e.
    // AlwaysUpdateMinDistance<true> always updates the given distance, while
    // AlwaysUpdateMinDistance<false> does not.  This optimization increases the
    // speed of GetDistance() by about 10% without creating code duplication.
    fun alwaysUpdateMinDistance(x: S2Point, a: S2Point, b: S2Point, min_dist: MutableS1ChordAngle, always_update: Boolean = true): Boolean {
        Assertions.assertPointIsUnitLength(x)
        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)

        val xa2 = (x - a).norm2()
        val xb2 = (x - b).norm2()
        if (alwaysUpdateMinInteriorDistance(x, a, b, xa2, xb2, min_dist, always_update)) {
            return true  // Minimum distance is attained along the edge interior.
        }
        // Otherwise the minimum distance is to one of the endpoints.
        val dist2 = min(xa2, xb2)
        if (!always_update && dist2 >= min_dist.length2) {
            return false
        }
        min_dist.length2 = dist2
        return true
    }

    // If the minimum distance from X to AB is attained at an interior point of AB
    // (i.e., not an endpoint), and that distance is less than "min_dist" or
    // "always_update" is true, then update "min_dist" and return true.  Otherwise
    // return false.
    //
    // The "Always" in the function name refers to the template argument, i.e.
    // AlwaysUpdateMinInteriorDistance<true> always updates the given distance,
    // while AlwaysUpdateMinInteriorDistance<false> does not.  This optimization
    // increases the speed of GetDistance() by about 10% without creating code
    // duplication.
    private fun alwaysUpdateMinInteriorDistance(x: S2Point, a: S2Point, b: S2Point, xa2: Double, xb2: Double, min_dist: MutableS1ChordAngle, always_update: Boolean = true): Boolean {
        Assertions.assertPointIsUnitLength(x)
        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)
        Assertions.assert { xa2 == (x - a).norm2() }
        Assertions.assert { xb2 == (x-b).norm2() }

        // The closest point on AB could either be one of the two vertices (the
        // "vertex case") or in the interior (the "interior case").  Let C = A x B.
        // If X is in the spherical wedge extending from A to B around the axis
        // through C, then we are in the interior case.  Otherwise we are in the
        // vertex case.
        //
        // Check whether we might be in the interior case.  For this to be true, XAB
        // and XBA must both be acute angles.  Checking this condition exactly is
        // expensive, so instead we consider the planar triangle ABX (which passes
        // through the sphere's interior).  The planar angles XAB and XBA are always
        // less than the corresponding spherical angles, so if we are in the
        // interior case then both of these angles must be acute.
        //
        // We check this by computing the squared edge lengths of the planar
        // triangle ABX, and testing acuteness using the law of cosines:
        //
        //             max(XA^2, XB^2) < min(XA^2, XB^2) + AB^2
        //
        if (max(xa2, xb2) >= min(xa2, xb2) + (a-b).norm2()) {
            return false
        }
        // The minimum distance might be to a point on the edge interior.  Let R
        // be closest point to X that lies on the great circle through AB.  Rather
        // than computing the geodesic distance along the surface of the sphere,
        // instead we compute the "chord length" through the sphere's interior.
        // If the squared chord length exceeds min_dist.length2() then we can
        // return "false" immediately.
        //
        // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
        // is the point X projected onto the plane through the great circle AB.
        // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
        // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
        // is faster and the corresponding distance on the Earth's surface is
        // accurate to within 1% for distances up to about 1800km.
        val c = S2.robustCrossProd(a, b)
        val c2 = c.norm2()
        val xDotC = x.dotProd(c)
        val xDotC2 = xDotC * xDotC
        if (!always_update && xDotC2 > c2 * min_dist.length2) {
            // The closest point on the great circle AB is too far away.  We need to
            // test this using ">" rather than ">=" because the actual minimum bound
            // on the distance is (x_dot_c2 / c2), which can be rounded differently
            // than the (more efficient) multiplicative test above.
            return false
        }
        // Otherwise we do the exact, more expensive test for the interior case.
        // This test is very likely to succeed because of the conservative planar
        // test we did initially.
        val cx = c.crossProd(x)
        if (a.dotProd(cx) >= 0 || b.dotProd(cx) <= 0) {
            return false
        }
        // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
        // This calculation has good accuracy for all chord lengths since it
        // is based on both the dot product and cross product (rather than
        // deriving one from the other).  However, note that the chord length
        // representation itself loses accuracy as the angle approaches Pi.
        val qr = 1 - sqrt(cx.norm2() / c2)
        val dist2 = (xDotC2 / c2) + (qr * qr)
        if (!always_update && dist2 >= min_dist.length2) {
            return false
        }
        min_dist.length2 = dist2
        return true
    }



}
