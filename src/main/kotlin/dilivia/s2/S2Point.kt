
package dilivia.s2

import com.google.common.geometry.S2LatLng
import kotlin.math.atan2

/**
 * An S2Point represents a point on the unit sphere as a 3D vector. Usually
 * points are normalized to be unit length, but some methods do not require
 * this.
 */
@Strictfp
class S2Point(coords: List<Double>) : RVector<S2Point, Double>(coords, DoubleType()) {

    @JvmOverloads
    constructor(x: Double = 0.0, y: Double = 0.0, z: Double= 0.0) : this(listOf<Double>(x, y, z)) {}

    constructor(x: Int, y: Int, z: Int) : this(x.toDouble(), y.toDouble(), z.toDouble())

    fun x(): Double {
        return get(0)
    }

    fun y(): Double {
        return get(1)
    }

    fun z(): Double {
        return get(2)
    }

    /**
     * return a vector orthogonal to this one
     */
    override fun ortho(): S2Point {
        val k = largestAbsComponent()
        val temp: S2Point
        temp = when (k) {
            1 -> S2Point(1, 0, 0)
            2 -> S2Point(0, 1, 0)
            else -> S2Point(0, 0, 1)
        }
        return normalize(crossProd(this, temp))
    }

    /**
     * Return the index of the largest component fabs
     */
    fun largestAbsComponent(): Int {
        val temp = abs(this)
        return if (temp.x() > temp.y()) {
            if (temp.x() > temp.z()) 0 else 2
        } else {
            if (temp.y() > temp.z()) 1 else 2
        }
    }

    /**
     * Return the angle between two vectors in radians
     */
    override fun angle(v: S2Point): Double {
        return atan2(crossProd(this, v).norm(), dotProd(v))
    }

    fun crossProd(other: S2Point): S2Point {
        return S2Point(
                y() * other.z() - z() * other.y(),
                z() * other.x() - x() * other.z(),
                x() * other.y() - y() * other.x()
        )
    }

    fun toDegreesString(): String {
        val s2LatLng = S2LatLng(this)
        return "(" + s2LatLng.latDegrees() + ", " + s2LatLng.lngDegrees() + ")"
    }

    /**
     * Calcualates hashcode based on stored coordinates. Since we want +0.0 and
     * -0.0 to be treated the same, we ignore the sign of the coordinates.
     *
     * @Override
     * public int hashCode() {
     * long value = 17;
     * value += 37 * value + Double.doubleToLongBits(Math.abs(x()));
     * value += 37 * value + Double.doubleToLongBits(Math.abs(y()));
     * value += 37 * value + Double.doubleToLongBits(Math.abs(z()));
     * return (int) (value ^ (value >>> 32));
     * }
     */
    override fun sqrt(): S2Point = S2Point(sqrtCoordinates())

    override fun normalize(): S2Point = S2Point(normalizedCoordinates())

    override operator fun plus(other: S2Point): S2Point = S2Point(sumCoordinates(other))

    override operator fun minus(other: S2Point): S2Point = S2Point(subtractCoordinates(other))

    override fun times(other: Double): S2Point = S2Point(scalarMulCoordinates(other))

    override fun div(other: Double): S2Point = S2Point(scalarDivCoordinates(other))

    override fun abs(): S2Point = S2Point(absCoordinates())

    companion object {

        @JvmStatic
        fun minus(p1: S2Point, p2: S2Point): S2Point = p1.minus(p2)

        @JvmStatic
        fun unaryMinus(p: S2Point): S2Point = p.unaryMinus()

        @JvmStatic
        fun crossProd(p1: S2Point, p2: S2Point): S2Point = p1.crossProd(p2)

        @JvmStatic
        fun plus(p1: S2Point, p2: S2Point): S2Point = p1.plus(p2)

        @JvmStatic
        fun times(p: S2Point, m: Double): S2Point = p.times(m)

        @JvmStatic
        fun div(p: S2Point, m: Double): S2Point = p.div(m)

        @JvmStatic
        fun abs(p: S2Point): S2Point = p.abs()

        @JvmStatic
        fun normalize(p: S2Point): S2Point {
            var norm = p.norm()
            if (norm != 0.0) {
                norm = 1.0 / norm
            }
            return times(p, norm)
        }
    }

}