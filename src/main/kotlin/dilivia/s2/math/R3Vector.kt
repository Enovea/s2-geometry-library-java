package dilivia.s2.math

import kotlin.math.atan2

@Strictfp
abstract class R3Vector<V, T> constructor(coords: List<T>, type: FloatingPointType<T>) : RVector<V, T>(coords, type) where V : R3Vector<V, T>, T: Number, T: Comparable<T> {

    init {
        require(coords.size == 3) { "Points must have exactly 3 coordinates" }
    }

    constructor(x: T, y: T, z: T, type: FloatingPointType<T>): this(listOf(x, y ,z), type)

    val x: T
            get() = coords[0]

    val y: T
        get()= coords[1]

    val z: T
            get() = coords[2]

    fun crossProd(other: V): V = newInstance(listOf(
            /*x =*/ type.minus(type.times(y,other.z), type.times(z, other.y)),
            /*y =*/ type.minus(type.times(z, other.x), type.times(x, other.z)),
            /*z =*/ type.minus(type.times(x, other.y), type.times(y, other.x))
    ))

    override fun ortho(): V {
        val k: Int = largestAbsComponent()
        val temp = when (k) {
            1 -> newInstance(listOf(type.one, type.zero, type.zero))
            2 -> newInstance(listOf(type.zero, type.one, type.zero))
            else -> newInstance(listOf(type.zero, type.zero, type.one))
        }
        return crossProd(temp).normalize()
    }

    override fun angle(v: V): T = type.atan2(crossProd(v).norm(), dotProd(v))

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
        fun <V, T> plus(p1: V, p2: V): V where V : R3Vector<V, T>, T: Number, T: Comparable<T> {
            return p1 + p2
        }

        @JvmStatic
        fun <V, T> times(p: V, m: T): V where V : R3Vector<V, T>, T: Number, T: Comparable<T> {
            return p * m
        }

        @JvmStatic
        fun <V, T> dotProd(p1: V, p2: V): T where V : R3Vector<V, T>, T: Number, T: Comparable<T> {
            return p1.dotProd(p2)
        }
    }
}