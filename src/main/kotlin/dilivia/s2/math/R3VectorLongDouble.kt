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

import ch.obermuhlner.math.big.BigDecimalMath
import java.math.BigDecimal
import java.math.MathContext
import java.math.RoundingMode


object LongDoubleType : FloatingPointType<BigDecimal> {

    val mathContext = MathContext(113, RoundingMode.HALF_EVEN)

    override val zero: BigDecimal = fromDouble(0.0)

    override val one: BigDecimal = fromDouble(1.0)

    override fun fromDouble(v: Double): BigDecimal {
        val bigDecimal = BigDecimal(v, mathContext)
        bigDecimal.setScale(2, RoundingMode.HALF_EVEN)
        return bigDecimal
    }

    override fun sqrt(v: BigDecimal): BigDecimal = v.sqrt(mathContext)

    override fun plus(a: BigDecimal, b: BigDecimal): BigDecimal = a.add(b, mathContext)

    override fun minus(a: BigDecimal, b: BigDecimal): BigDecimal = a.subtract(b, mathContext)

    override fun inv(a: BigDecimal): BigDecimal = one.divide(a, mathContext)

    override fun times(a: BigDecimal, b: BigDecimal): BigDecimal = a.multiply(b, mathContext)

    override fun div(a: BigDecimal, b: BigDecimal): BigDecimal = a.divide(b, mathContext)

    override fun abs(a: BigDecimal): BigDecimal = a.abs(mathContext)

    override fun sum(values: List<BigDecimal>): BigDecimal = values.reduce { acc, bigDecimal -> plus(acc, bigDecimal) }

    override fun equals(v1: BigDecimal, v2: BigDecimal): Boolean = v1 == v2

    override fun atan2(y: BigDecimal, x: BigDecimal): BigDecimal = BigDecimalMath.atan2(y, x, mathContext)
}

@Strictfp
open class R3VectorLongDouble constructor(coords: List<BigDecimal>) : R3Vector<R3VectorLongDouble, BigDecimal>(coords, LongDoubleType) {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    @JvmOverloads constructor(x: BigDecimal = LongDoubleType.zero, y: BigDecimal = LongDoubleType.zero, z: BigDecimal = LongDoubleType.zero): this(listOf(x, y ,z))

    @JvmOverloads constructor(x: Double, y: Double, z: Double): this(listOf(LongDoubleType.fromDouble(x), LongDoubleType.fromDouble(y), LongDoubleType.fromDouble(z)))

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(LongDoubleType.fromDouble(x.toDouble()), LongDoubleType.fromDouble(y.toDouble()), LongDoubleType.fromDouble(z.toDouble()))

    override fun newInstance(coords: List<BigDecimal>): R3VectorLongDouble = R3VectorLongDouble(coords)

    companion object {
        @JvmStatic
        fun plus(p1: R3VectorLongDouble, p2: R3VectorLongDouble): R3VectorLongDouble {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R3VectorLongDouble, m: BigDecimal): R3VectorLongDouble {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R3VectorLongDouble, p2: R3VectorLongDouble): BigDecimal {
            return p1.dotProd(p2)
        }
    }
}

fun Int.toLD(): BigDecimal = LongDoubleType.fromDouble(this.toDouble())
fun Double.toLD(): BigDecimal = LongDoubleType.fromDouble(this)