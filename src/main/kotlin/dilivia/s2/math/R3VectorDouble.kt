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

@Strictfp
open class R3VectorDouble constructor(coords: List<Double>) : R3Vector<R3VectorDouble, Double>(coords.map { if (it == -0.0) 0.0 else it }, DoubleType()) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0): this(listOf(x, y ,z))

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    override fun newInstance(coords: List<Double>): R3VectorDouble = R3VectorDouble(coords)

    companion object {
        @JvmStatic
        fun plus(p1: R3VectorDouble, p2: R3VectorDouble): R3VectorDouble {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R3VectorDouble, m: Double): R3VectorDouble {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R3VectorDouble, p2: R3VectorDouble): Double {
            return p1.dotProd(p2)
        }
    }
}