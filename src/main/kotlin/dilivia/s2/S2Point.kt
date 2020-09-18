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

import Matrix3x3
import com.google.common.geometry.S2
import dilivia.s2.Assertions.assertPointIsUnitLength
import dilivia.s2.math.*
import kotlin.math.atan2

/**
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually
 * points are normalized to be unit length, but some methods do not require
 * this.
 */
@Strictfp
open class S2Point(coords: List<Double>) : R3Vector<S2Point, Double>(coords.map { if (it == -0.0) 0.0 else it }, DoubleType()) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads
    constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0) : this(listOf<Double>(x, y, z)) {}

    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    constructor(vectorExactFloat: R3VectorExactFloat): this(vectorExactFloat.coords.map { it.toDouble() })

    fun x(): Double {
        return get(0)
    }

    fun y(): Double {
        return get(1)
    }

    fun z(): Double {
        return get(2)
    }

    fun toMutable(): MutableS2Point = MutableS2Point(coords.toMutableList())

    override fun newInstance(coords: List<Double>): S2Point = S2Point(coords)

    /**
     * return a vector orthogonal to this one
     */
    @Strictfp
    override fun ortho(): S2Point {
        val k = largestAbsComponent()
        val temp: S2Point
        temp = when (k) {
            1 -> S2Point(1, 0, 0)
            2 -> S2Point(0, 1, 0)
            else -> S2Point(0, 0, 1)
        }
        return normalize(crossProd(this, temp))
    }

    @Strictfp
    fun toDegreesString(): String {
        val s2LatLng = S2LatLng.fromPoint(this)
        return "(" + s2LatLng.latDegrees() + ", " + s2LatLng.lngDegrees() + ")"
    }

    /**
     * Return true if the given point is approximately unit length (this is mainly
     * useful for assertions).
     */
    @Strictfp
    fun isUnitLength(): Boolean = kotlin.math.abs(norm2() - 1) <= 5 * S2.DBL_EPSILON

    @Strictfp
    fun toLongDouble(): R3VectorLongDouble = R3VectorLongDouble(x, y, z)

    @Strictfp
    fun toExactFloat(): R3VectorExactFloat = R3VectorExactFloat(x, y, z)

    companion object {

        /**
         * Return a unit-length vector that is orthogonal to "a". Satisfies Ortho(-a)
         * = -Ortho(a) for all a.
         */
        @JvmStatic
        @Strictfp
        fun ortho(a: S2Point): S2Point {
            // The current implementation in S2Point has the property we need,
            // i.e. Ortho(-a) = -Ortho(a) for all a.
            return a.ortho()
        }

        /**
         * Return a unique "origin" on the sphere for operations that need a fixed
         * reference point. It should *not* be a point that is commonly used in edge
         * tests in order to avoid triggering code to handle degenerate cases. (This
         * rules out the north and south poles.)
         */
        @JvmStatic
        @Strictfp
        fun origin(): S2Point {
            return S2Point(0, 1, 0)
        }

        @JvmStatic
        @Strictfp
        fun isUnitLength(p: S2Point): Boolean = p.isUnitLength()

        @JvmStatic
        @Strictfp
        fun minus(p1: S2Point, p2: S2Point): S2Point = p1.minus(p2)

        @JvmStatic
        @Strictfp
        fun unaryMinus(p: S2Point): S2Point = p.unaryMinus()

        @JvmStatic
        @Strictfp
        fun crossProd(p1: S2Point, p2: S2Point): S2Point = p1.crossProd(p2)

        @JvmStatic
        @Strictfp
        fun plus(p1: S2Point, p2: S2Point): S2Point = p1.plus(p2)

        @JvmStatic
        @Strictfp
        fun times(p: S2Point, m: Double): S2Point = p.times(m)

        @JvmStatic
        @Strictfp
        fun div(p: S2Point, m: Double): S2Point = p.div(m)

        @JvmStatic
        @Strictfp
        fun abs(p: S2Point): S2Point = p.abs()

        @JvmStatic
        @Strictfp
        fun normalize(p: S2Point): S2Point {
            var norm = p.norm()
            if (norm != 0.0) {
                norm = 1.0 / norm
            }
            return times(p, norm)
        }

        /**
         * Return true if two points are within the given distance of each other
         * (mainly useful for testing).
         */
        //@JvmStatic
        //fun approxEquals(a: S2Point, b: S2Point, maxError: Double): Boolean = a.angle(b) <= maxError

        //@JvmStatic
        //fun approxEquals(a: S2Point, b: S2Point): Boolean = approxEquals(a, b, 1e-15)

        @JvmStatic
        @JvmOverloads
        @Strictfp
        fun approxEquals(a: S2Point, b: S2Point, maxErrorAngle: S1Angle= S1Angle.radians(1e-15)): Boolean = S1Angle(a, b) <= maxErrorAngle

        @JvmStatic
        @Strictfp
        fun getFrame(z: S2Point): Matrix3x3 {
            assertPointIsUnitLength(z)
            val m = Matrix3x3()
            m.setCol(2, z)
            m.setCol(1, ortho(z))
            m.setCol(0, m.col(1).crossProd(z));  // Already unit-length.
            return m
        }

        @JvmStatic
        @Strictfp
        fun toFrame(m: Matrix3x3, p: S2Point): S2Point {
            // The inverse of an orthonormal matrix is its transpose.
            return m.transpose() * p;
        }

        @JvmStatic
        @Strictfp
        fun fromFrame(m: Matrix3x3, q: S2Point): S2Point {
            return m * q;
        }

        /**
        // Return true if the points A, B, C are strictly counterclockwise.  Return
        // false if the points are clockwise or collinear (i.e. if they are all
        // contained on some great circle).
        //
        // Due to numerical errors, situations may arise that are mathematically
        // impossible, e.g. ABC may be considered strictly CCW while BCA is not.
        // However, the implementation guarantees the following:
        //
        //   If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
                ABSL_DEPRECATED("Use s2pred::Sign instead.")
         */
        @Deprecated("Use S2Predicates.sign instead.", replaceWith = ReplaceWith("S2Predicates.sign"))
        @JvmStatic
        @Strictfp
        fun simpleCCW(a: S2Point, b: S2Point, c: S2Point): Boolean {
            // We compute the signed volume of the parallelepiped ABC.  The usual
            // formula for this is (AxB).C, but we compute it here using (CxA).B
            // in order to ensure that ABC and CBA are not both CCW.  This follows
            // from the following identities (which are true numerically, not just
            // mathematically):
            //
            //     (1) x.CrossProd(y) == -(y.CrossProd(x))
            //     (2) (-x).DotProd(y) == -(x.DotProd(y))

            return c.crossProd(a).dotProd(b) > 0;
        }

        // Return a vector "c" that is orthogonal to the given unit-length vectors
        // "a" and "b".  This function is similar to a.CrossProd(b) except that it
        // does a better job of ensuring orthogonality when "a" is nearly parallel
        // to "b", and it returns a non-zero result even when a == b or a == -b.
        //
        // It satisfies the following properties (RCP == RobustCrossProd):
        //
        //   (1) RCP(a,b) != 0 for all a, b
        //   (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b
        //   (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
        //   (4) RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
        //
        // The result is not guaranteed to be unit length.
        @JvmStatic
        @Strictfp
        fun robustCrossProd(a: S2Point, b: S2Point): S2Point {
            // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
            // approaches zero.  This leads to situations where a.CrossProd(b) is not
            // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
            // but we also want b.RobustCrossProd(a) == -a.RobustCrossProd(b).
            //
            // The easiest fix is to just compute the cross product of (b+a) and (b-a).
            // Mathematically, this cross product is exactly twice the cross product of
            // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
            // are always perpendicular (since "a" and "b" are unit length).  This
            // yields a result that is nearly orthogonal to both "a" and "b" even if
            // these two values differ only in the lowest bit of one component.
            assertPointIsUnitLength(a)
            assertPointIsUnitLength(b)
            val x = (b + a).crossProd(b - a);
            if (x != S2Point(0, 0, 0)) return x

            // The only result that makes sense mathematically is to return zero, but
            // we find it more convenient to return an arbitrary orthogonal vector.
            return ortho(a)
        }

    }

}

operator fun Double.times(point: S2Point): S2Point = point * this