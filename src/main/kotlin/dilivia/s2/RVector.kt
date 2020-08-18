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

import kotlin.math.abs
import kotlin.math.atan2
import kotlin.math.sqrt

interface FloatingPointType<T: Number> {

    fun sqrt(v: T): T

    fun plus(a: T, b: T): T

    fun minus(a: T, b: T): T

    fun inv(a: T): T

    fun times(a: T, b: T): T

    fun div(a: T, b: T): T

    fun abs(a: T): T

    fun sum(values: List<T>): T
}

class DoubleType() : FloatingPointType<Double> {
    override fun sqrt(v: Double): Double = sqrt(v)

    override fun plus(a: Double, b: Double): Double = a + b

    override fun minus(a: Double, b: Double): Double = a - b

    override fun inv(a: Double): Double = 1.0 / a

    override fun times(a: Double, b: Double): Double = a * b

    override fun div(a: Double, b: Double): Double = a / b

    override fun abs(a: Double): Double = abs(a)

    override fun sum(values: List<Double>): Double = values.sum()

}

/**
 *
 */
abstract class RVector<V, T>(val coords: List<T>, val type: FloatingPointType<T>) : Comparable<RVector<V, T>> where V : RVector<V, T>, T : Number,  T: Comparable<T> {

    val size: Int
        get() = coords.size


    operator fun get(index: Int): T {
        check(index in 0 until size) { "Coordinate index $index is not in range 0..${size - 1}" }
        return coords[index]
    }

    fun norm2(): T {
        return this.dotProd(this)
    }

    fun norm(): T = type.sqrt(norm2())

    abstract fun angle(v: V): Double

    abstract fun ortho(): V

    abstract fun sqrt(): V

    protected fun sqrtCoordinates() = coords.map { type.sqrt(it) }

    // Normalized vector if the norm is nonzero. Not for integer types.
    abstract fun normalize(): V

    protected fun normalizedCoordinates(): List<T> {
        var n = norm()
        check(n != 0.0)
        n = type.inv(n)
        return coords.map { type.times(it, n) }
    }

    override operator fun compareTo(other: RVector<V, T>): Int {
        require(size == other.size)
        for (i in 0 until size) {
            val compareTo = this[i].compareTo(other[i])
            if (i == coords.lastIndex || compareTo != 0) return compareTo
        }
        throw IllegalStateException("Return should have been called in the for loop.")
    }

    abstract operator fun plus(other: V): V

    protected fun sumCoordinates(other: V): List<T> {
        require(size == other.size)
        return coords.indices.map { i -> type.plus(this[i], other[i]) }
    }

    abstract operator fun minus(other: V): V

    protected fun subtractCoordinates(other: V): List<T> {
        require(size == other.size)
        return coords.indices.map { i -> type.minus(this[i], other[i]) }
    }

    abstract operator fun times(other: Double): V

    protected fun scalarMulCoordinates(other: T): List<T> = coords.indices.map { i -> type.times(this[i], other) }

    abstract operator fun div(other: Double): V

    protected fun scalarDivCoordinates(other: T) = coords.indices.map { i -> type.div(this[i], other) }

    fun dotProd(other: RVector<V, T>): T {
        return type.sum(coords.indices.map { type.times(this[it], other[it]) })
    }

    abstract fun abs(): V

    protected fun absCoordinates(): List<T> {
        return coords.indices.map { i -> type.abs(this[i]) }
    }

    override fun equals(other: Any?): Boolean {
        if (!this::class.isInstance(other)) {
            return false
        }
        return this.size == (other as V).size && this.coords.indices.all { i -> this[i] == other[i] }
    }

    fun approxEquals(other: V, margin: Double = 1e-15): Boolean {
        return coords.indices.all { i -> type.abs(type.minus(this[i], other[i])).toDouble() < margin }
    }

    /**
     * Calculates hashcode based on stored coordinates. Since we want +0.0 and
     * -0.0 to be treated the same, we ignore the sign of the coordinates.
     */
    override fun hashCode(): Int {
        var value: Long = 17
        coords.forEach { c ->
            value += 37 * value + java.lang.Double.doubleToLongBits(type.abs(c).toDouble())
        }
        return (value xor (value ushr 32)).toInt()
    }

    override fun toString(): String {
        return coords.joinToString(prefix = "(", postfix = ")")
    }
}

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
class R2Vector @JvmOverloads constructor(x: Double = 0.0, y: Double = 0.0) : RVector<R2Vector, Double>(listOf(x, y), DoubleType()) {

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

    fun x(): Double = coords[0]

    fun y(): Double = coords[1]

    override fun sqrt(): R2Vector = R2Vector(sqrtCoordinates())

    override fun normalize(): R2Vector = R2Vector(normalizedCoordinates())

    override operator fun plus(other: R2Vector): R2Vector = R2Vector(sumCoordinates(other))

    override operator fun minus(other: R2Vector): R2Vector = R2Vector(subtractCoordinates(other))

    override operator fun times(other: Double): R2Vector = R2Vector(scalarMulCoordinates(other))

    override operator fun div(other: Double): R2Vector = R2Vector(scalarDivCoordinates(other))

    fun crossProd(other: R2Vector): Double = x() * other.y() - y() * other.x()

    // Returns the angle between "this" and v in radians. If either vector is
    // zero-length, or nearly zero-length, the result will be zero, regardless of
    // the other value.
    override fun angle(v: R2Vector): Double = atan2(crossProd(v), dotProd(v))

    // return a vector orthogonal to the current one
    // with the same norm and counterclockwise to it
    override fun ortho(): R2Vector = R2Vector(-y(), x())

    override fun abs(): R2Vector = R2Vector(absCoordinates())

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

@Strictfp
class R3Vector @JvmOverloads constructor(x: Double = 0.0, y: Double = 0.0, z: Double = 0.0) : RVector<R3Vector, Double>(listOf(x, y, z), DoubleType()) {

    /**
     * Create a R3Vector instance from Int values.
     *
     * @param x x coordinate as a integer.
     * @param y y coordinate as a integer.
     * @param z z coordinate as a integer.
     */
    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    constructor(coord: List<Double>) : this(coord[0], coord[1], coord[2]) {
        require(coord.size == 3) { "Points must have exactly 3 coordinates" }
    }

    fun x(): Double = coords[0]

    fun y(): Double = coords[1]

    fun z(): Double = coords[2]

    override fun sqrt(): R3Vector = R3Vector(sqrtCoordinates())

    override fun normalize(): R3Vector = R3Vector(normalizedCoordinates())

    override operator fun plus(other: R3Vector): R3Vector = R3Vector(sumCoordinates(other))

    override operator fun minus(other: R3Vector): R3Vector = R3Vector(subtractCoordinates(other))

    override operator fun times(other: Double): R3Vector = R3Vector(scalarMulCoordinates(other))

    override operator fun div(other: Double): R3Vector = R3Vector(scalarDivCoordinates(other))

    fun crossProd(other: R3Vector): R3Vector = R3Vector(
            x = y() * other.z() - z() * other.y(),
            y = z() * other.x() - x() * other.z(),
            z = x() * other.y() - y() * other.x()
    )

    override fun ortho(): R3Vector {
        val k: Int = largestAbsComponent()
        val temp = when (k) {
            1 -> R3Vector(1, 0, 0)
            2 -> R3Vector(0, 1, 0)
            else -> R3Vector(0, 0, 1)
        }
        return crossProd(temp).normalize()
    }

    override fun angle(v: R3Vector): Double = atan2(crossProd(v).norm(), dotProd(v))

    override fun abs(): R3Vector = R3Vector(absCoordinates())

    // return the index of the largest component (fabs)
    fun largestAbsComponent(): Int {
        val temp = abs()
        return if (temp[0] > temp[1]) {
            if (temp[0] > temp[2]) 0 else 2
        } else {
            if (temp[1] > temp[2]) 1 else 2
        }
    }

    companion object {
        @JvmStatic
        fun plus(p1: R3Vector, p2: R3Vector): R3Vector {
            return p1 + p2
        }

        @JvmStatic
        fun times(p: R3Vector, m: Double): R3Vector {
            return p * m
        }

        @JvmStatic
        fun dotProd(p1: R3Vector, p2: R3Vector): Double {
            return p1.dotProd(p2)
        }
    }
}



