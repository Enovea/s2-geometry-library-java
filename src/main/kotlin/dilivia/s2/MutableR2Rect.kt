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

import dilivia.s2.math.R2Point

@Strictfp
class MutableR2Rect(x: MutableR1Interval, y: MutableR1Interval) : R2Rect(x, y) {

    constructor(): this(MutableR1Interval.empty(), MutableR1Interval.empty())

    override operator fun get(idx: Int): MutableR1Interval {
        require(idx in 0..1)
        return if (idx == 0) x as MutableR1Interval else y as MutableR1Interval
    }

    operator fun set(idx: Int, v: R1Interval) {
        require(idx in 0..1)
        if (idx == 0) {
            (x as MutableR1Interval)[0] = v[0]
            x[1] = v[1]
        } else {
            (y as MutableR1Interval)[0] = v[0]
            y[1] = v[1]
        }
    }

    override fun addPoint(p: R2Point): MutableR2Rect {
        x.addPoint(p[0])
        y.addPoint(p[1])
        return this
    }

    // Expand the rectangle to include the given other rectangle.  This is the
    // same as replacing the rectangle by the union of the two rectangles, but
    // is somewhat more efficient.
    override fun addRect(other: R2Rect): MutableR2Rect {
        (x as MutableR1Interval).addInterval(other.x)
        (y as MutableR1Interval).addInterval(other.y)
        return this
    }

    fun expands(margin: R2Point): MutableR2Rect {
        (x as MutableR1Interval).expands(margin.x())
        (y as MutableR1Interval).expands(margin.y())
        if (x.isEmpty || y.isEmpty) {
            x.setEmpty()
            y.setEmpty()
        }
        return this
    }

    fun expands(margin: Double): MutableR2Rect = expands(R2Point(margin, margin))

    fun setUnion(other: R2Rect): MutableR2Rect {
        (x as MutableR1Interval).setUnion(other.x)
        (y as MutableR1Interval).setUnion(other.y)
        return this
    }

    fun setIntersection(other: R2Rect): MutableR2Rect {
        (x as MutableR1Interval).setIntersection(other.x)
        (y as MutableR1Interval).setIntersection(other.y)
        if (x.isEmpty || y.isEmpty) {
            x.setEmpty()
            y.setEmpty()
        }
        return this
    }
}