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

/**
 * An R2Rect represents a closed axis-aligned rectangle in the (x,y) plane.
 *
 * @property x Interval in x.
 * @property y Interval in y.
 * @constructor Construct a rectangle from the given intervals in x and y. The two intervals must either be both empty
 * or both non-empty.
 * @param x Interval in x.
 * @param y Interval in y.
 *
 * @author Fabien Meurisse <fabien.meurisse@enovea.net>
 * @since 1.0
 */
open class R2Rect(open val x: R1Interval, open val y: R1Interval) {

    init {
        check(isValid)
    }

    /**
     * Construct a rectangle from the given lower-left and upper-right points.
     *
     * @param lo Lower-left corner of the rectangle.
     * @param hi Upper-Right corner of the rectangle.
     */
    constructor(lo: R2Point, hi: R2Point) : this(R1Interval(lo.x(), hi.x()), R1Interval(lo.y(), hi.y()))

    /**
     * The default constructor creates an empty R2Rect.
     */
    constructor() : this(R1Interval.empty(), R1Interval.empty())

    /** Lower-left corner. */
    val lo: R2Point
        get() = R2Point(x.lo, y.lo)

    /** Upper-right corner. */
    val hi: R2Point
        get() = R2Point(x.hi, y.hi)

    /**
     * Get the i-th interval (x = 0, y = 1) of the rectangle.
     *
     * @param idx The interval index.
     * @return x for idx = 0 and y for idx = 1
     */
    operator fun get(idx: Int): R1Interval {
        require(idx in 0..1)
        return if (idx == 0) x else y
    }

    /**
     * Indicates if the rectangle is valid, which essentially just means that if the bound for either axis is empty
     * then both must be.
     */
    val isValid: Boolean
        get() {
            return x.isEmpty == y.isEmpty // The x/y ranges must either be both empty or both non-empty.
        }

    /** Indicates if the rectangle is empty, i.e. it contains no points at all. */
    val isEmpty: Boolean
        get() = x.isEmpty

    /**
     * Get the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.
     * Vertex 0 is in the lower-left corner.  For convenience, the argument is reduced modulo 4 to the range [0..3].
     *
     * @param k The counter-clockwise index of the vertex.
     * @return The k-th vertex.
     */
    fun getVertex(k: Int): R2Point {
        // Twiddle bits to return the points in CCW order (lower left, lower right, upper right, upper left).
        val j = (k shr 1) and 1
        return getVertex(j xor (k and 1), j)
    }

    /**
    // Return the vertex in direction "i" along the x-axis (0=left, 1=right) and
    // direction "j" along the y-axis (0=down, 1=up).  Equivalently, return the
    // vertex constructed by selecting endpoint "i" of the x-interval (0=lo,
    // 1=hi) and vertex "j" of the y-interval.
    */
    fun getVertex(i: Int, j: Int): R2Point = R2Point(this[0][i], this[1][j])

    /** the center of the rectangle in (x,y)-space. */
    val center: R2Point
        get() = R2Point(x.center, y.center)

    /** the width and height of this rectangle in (x,y)-space. Empty rectangles have a negative width and height. */
    val size: R2Vector
        get() = R2Vector(x.length, y.length)

    // Return true if the rectangle contains the given point.  Note that
    // rectangles are closed regions, i.e. they contain their boundary.
    fun contains(p: R2Point): Boolean = x.contains(p.x()) && y.contains(p.y())

    // Return true if and only if the given point is contained in the interior
    // of the region (i.e. the region excluding its boundary).
    fun interiorContains(p: R2Point): Boolean = x.interiorContains(p.x()) && y.interiorContains(p.y())

    // Return true if and only if the rectangle contains the given other
    // rectangle.
    fun contains(other: R2Rect): Boolean = x.contains(other.x) && y.contains(other.y)

    // Return true if and only if the interior of this rectangle contains all
    // points of the given other rectangle (including its boundary).
    fun interiorContains(other: R2Rect): Boolean = x.interiorContains(other.x) && y.interiorContains(other.y)

    // Return true if this rectangle and the given other rectangle have any
    // points in common.
    fun intersects(other: R2Rect): Boolean = x.intersects(other.x) && y.intersects(other.y)

    // Return true if and only if the interior of this rectangle intersects
    // any point (including the boundary) of the given other rectangle.
    fun interiorIntersects(other: R2Rect): Boolean = x.interiorIntersects(other.x) && y.interiorIntersects(other.y)

    // Expand the rectangle to include the given point.  The rectangle is
    // expanded by the minimum amount possible.
    open fun addPoint(p: R2Point): R2Rect = R2Rect(
            x.addPoint(p[0]),
            y.addPoint(p[1])
    )

    // Expand the rectangle to include the given other rectangle.  This is the
    // same as replacing the rectangle by the union of the two rectangles, but
    // is somewhat more efficient.
    open fun addRect(other: R2Rect) = R2Rect(
            x.addInterval(other.x),
            y.addInterval(other.y)
    )

    // Return the closest point in the rectangle to the given point "p".
    // The rectangle must be non-empty.
    fun project(p: R2Point): R2Point = R2Point(x.project(p.x()), y.project(p.y()))

    // Return a rectangle that has been expanded on each side in the x-direction
    // by margin.x(), and on each side in the y-direction by margin.y().  If
    // either margin is empty, then shrink the interval on the corresponding
    // sides instead.  The resulting rectangle may be empty.  Any expansion of
    // an empty rectangle remains empty.
    fun expanded(margin: R2Point): R2Rect {
        val xx = x.expanded(margin.x())
        val yy = y.expanded(margin.y())
        if (xx.isEmpty || yy.isEmpty) return empty()
        return R2Rect(xx, yy)
    }

    fun expanded(margin: Double): R2Rect = expanded(R2Point(margin, margin))

    // Return the smallest rectangle containing the union of this rectangle and
    // the given rectangle.
    fun union(other: R2Rect): R2Rect = R2Rect(x.union(other.x), y.union(other.y))

    // Return the smallest rectangle containing the intersection of this
    // rectangle and the given rectangle.
    fun intersection(other: R2Rect): R2Rect {
        val xx = x.intersection(other.x)
        val yy = y.intersection(other.y)
        if (xx.isEmpty || yy.isEmpty) return empty()
        return R2Rect(xx, yy)
    }

    override fun equals(other: Any?): Boolean {
        if (other !is R2Rect) return false
        return x == other.x && y == other.y
    }

    override fun hashCode(): Int {
        var result = x.hashCode()
        result = 31 * result + y.hashCode()
        return result
    }

    // Return true if the x- and y-intervals of the two rectangles are the same
    // up to the given tolerance (see r1interval.h for details).
    fun approxEquals(other: R2Rect, maxError: Double = 1e-15): Boolean = (x.approxEquals(other.x, maxError) && y.approxEquals(other.y, maxError))

    override fun toString(): String = "[Lo $lo, Hi $hi]"

    companion object {

        /**
         * The canonical empty rectangle. Use isEmpty() to test for empty rectangles, since they have more than one
         * representation.
         */
        @JvmStatic
        fun empty() = R2Rect()

        /**
         * Construct a rectangle from a center point and size in each dimension.
         * Both components of size should be non-negative, i.e. this method cannot be used to create an empty rectangle.
         */
        @JvmStatic
        fun fromCenterSize(center: R2Point, size: R2Point): R2Rect = R2Rect(
                R1Interval(center.x() - 0.5 * size.x(), center.x() + 0.5 * size.x()),
                R1Interval(center.y() - 0.5 * size.y(), center.y() + 0.5 * size.y())
        )

        /**
         * Convenience method to construct a rectangle containing a single point.
         */
        @JvmStatic
        fun fromPoint(p: R2Point) = R2Rect(p, p)

        /**
         * Convenience method to construct the minimal bounding rectangle containing the two given points.
         * This is equivalent to starting with an empty rectangle and calling AddPoint() twice.  Note that it is
         * different than the R2Rect(lo, hi) constructor, where the first point is always used as the lower-left corner
         * of the resulting rectangle.
         */
        @JvmStatic
        fun fromPointPair(p1: R2Point, p2: R2Point): R2Rect = R2Rect(
                R1Interval.fromPointPair(p1.x(), p2.x()),
                R1Interval.fromPointPair(p1.y(), p2.y())
        )

    }
}

