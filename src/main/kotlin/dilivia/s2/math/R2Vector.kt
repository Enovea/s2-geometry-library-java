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
package dilivia.s2.math

import dilivia.s2.MutableR2Point
import kotlin.math.atan2

/**
 * R2Vector represents a vector in the two-dimensional space. It defines the basic geometrical operations for 2D
 * vectors, e.g. cross product, addition, norm, comparison etc.
 *
 * @property x First coordinate of the vector.
 * @property y Second coordinate of the vector.
 * @constructor Create a R2Vector instance.
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
@Strictfp
open class R2Vector @JvmOverloads constructor(x: Double = 0.0, y: Double = 0.0) : RVector<R2Vector, Double>(listOf(x, y), DoubleType()) {

    /**
     * Create a R2Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     */
    constructor(x: Int, y: Int) : this(x.toDouble(), y.toDouble())

    constructor(coord: List<Double>) : this(coord[0], coord[1]) {
        require(coord.size == 2) { "Points must have exactly 2 coordinates" }
    }

    constructor(coord: DoubleArray): this(coord[0], coord[1]) {
        require(coord.size == 2) { "Points must have exactly 2 coordinates" }
    }

    fun x(): Double = coords[0]

    fun y(): Double = coords[1]

    override fun newInstance(coords: List<Double>): R2Vector = R2Vector(coords)

    fun crossProd(other: R2Vector): Double = x() * other.y() - y() * other.x()

    // Returns the angle between "this" and v in radians. If either vector is
    // zero-length, or nearly zero-length, the result will be zero, regardless of
    // the other value.
    override fun angle(v: R2Vector): Double = atan2(crossProd(v), dotProd(v))

    // return a vector orthogonal to the current one
    // with the same norm and counterclockwise to it
    override fun ortho(): R2Vector = R2Vector(-y(), x())

    fun toMutable(): MutableR2Point = MutableR2Point(x(), y())

    companion object {
        @JvmStatic
        fun plus(p1: R2Vector, p2: R2Vector): R2Vector {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R2Vector, m: Double): R2Vector {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R2Vector, p2: R2Vector): Double {
            return p1.dotProd(p2)
        }
    }
}

typealias R2Point = R2Vector

operator fun Double.times(p: R2Vector): R2Vector = p * this