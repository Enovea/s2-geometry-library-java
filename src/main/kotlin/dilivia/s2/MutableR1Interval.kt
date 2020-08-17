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

import kotlin.math.max
import kotlin.math.min

/**
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
class MutableR1Interval(override var lo: Double, override var hi: Double): R1Interval(lo, hi) {

    constructor(lo: Int, hi: Int): this(lo.toDouble(), hi.toDouble())

    operator fun set(i: Int, v: Double) {
        when(i) {
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

    override fun addPoint(p: Double): R1Interval {
        when {
            isEmpty -> {
                lo = p
                hi = p
            }
            p < lo -> lo = p
            p > hi -> hi = p
        }
        return this
    }

    override fun addInterval(y: R1Interval): R1Interval {
        if (y.isNotEmpty) {
            if (isEmpty) {
                lo = y.lo
                hi = y.hi
            }
            else {
                lo = min(lo, y.lo)
                hi = max(hi, y.hi)
            }
        }
        return this
    }

    fun expands(radius: Double): R1Interval {
        // assert (radius >= 0);
        if (isNotEmpty) {
            lo -= radius
            hi += radius
        }
        return this
    }

    fun setUnion(y: R1Interval): R1Interval {
        if (isEmpty) {
            lo = y.lo
            hi = y.hi
        }
        else if (y.isNotEmpty) {
            lo = min(lo, y.lo)
            hi = max(hi, y.hi)
        }
        return this
    }

    fun setIntersection(y: R1Interval): R1Interval {
        lo = max(lo, y.lo)
        hi = min(hi, y.hi)
        return this
    }

}