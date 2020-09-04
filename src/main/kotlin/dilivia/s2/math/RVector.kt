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

/**
 *
 */
@Strictfp
abstract class RVector<V, T>(open val coords: List<T>, val type: FloatingPointType<T>) : Comparable<RVector<V, T>> where V : RVector<V, T>, T : Number, T: Comparable<T> {

    val size: Int
        get() = coords.size

    protected abstract fun newInstance(coords: List<T>): V

    operator fun get(index: Int): T {
        check(index in 0 until size) { "Coordinate index $index is not in range 0..${size - 1}" }
        return coords[index]
    }

    fun norm2(): T {
        return this.dotProd(this)
    }

    fun norm(): T = type.sqrt(norm2())

    abstract fun angle(v: V): T

    abstract fun ortho(): V

    fun sqrt(): V = newInstance(sqrtCoordinates())

    protected fun sqrtCoordinates() = coords.map { type.sqrt(it) }

    // Normalized vector if the norm is nonzero. Not for integer types.
    fun normalize(): V = newInstance(normalizedCoordinates())

    protected fun normalizedCoordinates(): List<T> {
        var n = norm()
        check(n != 0.0) { "|$this| = 0" }
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

    operator fun plus(other: V): V = newInstance(sumCoordinates(other))

    protected fun sumCoordinates(other: V): List<T> {
        require(size == other.size)
        return coords.indices.map { i -> type.plus(this[i], other[i]) }
    }

    operator fun minus(other: V): V = newInstance(subtractCoordinates(other))

    protected fun subtractCoordinates(other: V): List<T> {
        require(size == other.size)
        return coords.indices.map { i -> type.minus(this[i], other[i]) }
    }

    operator fun times(other: T): V = newInstance(scalarMulCoordinates(other))

    operator fun unaryMinus(): V = this * type.fromDouble(-1.0)

    operator fun times(other: V): V = newInstance(mulComponents(other))

    protected fun mulComponents(other: RVector<V, T>): List<T> = coords.indices.map { i -> type.times(this[i], other[i]) }

    protected fun scalarMulCoordinates(other: T): List<T> = coords.indices.map { i -> type.times(this[i], other) }

    operator fun div(other: T): V = newInstance(scalarDivCoordinates(other))

    protected fun scalarDivCoordinates(other: T) = coords.indices.map { i -> type.div(this[i], other) }

    fun dotProd(other: RVector<V, T>): T {
        return type.sum(coords.indices.map { type.times(this[it], other[it]) })
    }

    fun abs(): V = newInstance(absCoordinates())

    protected fun absCoordinates(): List<T> {
        return coords.indices.map { i -> type.abs(this[i]) }
    }

    override fun equals(other: Any?): Boolean {
        if (!this::class.isInstance(other)) {
            return false
        }
        return this.size == (other as V).size && this.coords.indices.all { i -> type.equals(this[i], other[i]) }
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




