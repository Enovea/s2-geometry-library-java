package dilivia.s2.math

interface FloatingPointType<T: Number> {

    val zero: T

    val one: T

    fun fromDouble(v: Double): T

    fun sqrt(v: T): T

    fun plus(a: T, b: T): T

    fun minus(a: T, b: T): T

    fun inv(a: T): T

    fun times(a: T, b: T): T

    fun div(a: T, b: T): T

    fun abs(a: T): T

    fun sum(values: List<T>): T

    fun equals(v1: T, v2: T): Boolean

    fun atan2(y: T, x: T): T

}