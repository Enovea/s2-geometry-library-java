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

import com.google.common.geometry.S2.M_PI
import dilivia.s2.math.R2Vector
import kotlin.math.abs

/**
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
@Strictfp
class MutableS1Interval(override var lo: Double, override var hi: Double) : S1Interval(lo, hi) {

    constructor(lo: Int, hi: Int) : this(lo.toDouble(), hi.toDouble())

    constructor(interval: S1Interval) : this(interval.lo, interval.hi)

    operator fun set(i: Int, v: Double) {
        when (i) {
            0 -> lo = v
            1 -> hi = v
            else -> throw ArrayIndexOutOfBoundsException("Index $i is out of bounds 0..1")
        }
    }

    override var bounds: R2Vector
        get() = super.bounds
        set(value) {
            lo = value[0]
            hi = value[1]
        }

    override fun addPoint(p: Double): MutableS1Interval {
        require(abs(p) <= M_PI)
        var cp = p
        if (cp == -M_PI) {
            cp = M_PI
        }
        if (!fastContains(cp)) {
            if (isEmpty) {
                lo = p
                hi = p
            } else {
                // Compute distance from p to each endpoint.
                val dlo = positiveDistance(cp, lo)
                val dhi = positiveDistance(hi, cp)
                if (dlo < dhi) {
                    lo = cp
                } else {
                    hi = cp
                }
                // Adding a point can never turn a non-full interval into a full one.
            }
        }
        return this
    }

    fun expands(radius: Double): MutableS1Interval {
        require(radius >= 0)
        if (isNotEmpty) {
            lo -= radius
            hi += radius
        }
        return this
    }

    fun setUnion(y: S1Interval): MutableS1Interval {
        TODO()
    }

    fun setIntersection(y: R1Interval): MutableS1Interval {
        TODO()
    }

    fun setEmpty(): MutableS1Interval {
        lo = M_PI
        hi = - M_PI
        return this
    }


    companion object {

        /**
         * Returns an empty interval. (Any interval where lo > hi is considered empty.)
         *
         * @return An empty interval instance with lo = 1 and hi = 0
         */
        @JvmStatic
        fun empty(): MutableS1Interval {
            return MutableS1Interval(empty)
        }

        /**
         * Convenience method to construct an interval containing a single point.
         *
         * @param p The unique point in the interval
         * @return An interval instance that contains only the given point.
         */
        @JvmStatic
        fun fromPoint(p: Double): MutableS1Interval {
            return MutableS1Interval(p, p)
        }

        /**
         * Convenience method to construct the minimal interval containing the two given points. This is equivalent to
         * starting with an empty interval and calling addPoint() twice, but it is more efficient.
         *
         * @param p1 A point to be included.
         * @param p2 A point to be included.
         * @return The minimum interval that contains the two points.
         */
        @JvmStatic
        fun fromPointPair(p1: Double, p2: Double): MutableS1Interval = TODO()

    }

}