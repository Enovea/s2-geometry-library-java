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

import com.google.common.geometry.S2
import com.google.common.geometry.S2Projections
import dilivia.s2.S1Interval.Companion.full
import dilivia.s2.S2Cap.Companion.fromCenterHeight
import dilivia.s2.S2CellId.Companion.fromLatLng
import dilivia.s2.S2CellId.Companion.fromPoint
import dilivia.s2.S2CellId.Companion.none
import dilivia.s2.S2LatLngRect.Companion.fullLat
import dilivia.s2.S2Point.Companion.crossProd
import dilivia.s2.S2Point.Companion.minus
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.S2Point.Companion.unaryMinus
import dilivia.s2.math.R2Vector

/**
 * An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it
 * supports efficient containment and intersection tests. However, it is also a
 * more expensive representation.
 */
@Strictfp
class S2Cell : S2Region {
    var face: Byte = 0
    var level: Byte = 0
    var orientation: Byte = 0
    var cellId: S2CellId = none()
    var uv = Array(2) { DoubleArray(2) }

    /**
     * Default constructor used only internally.
     */
    constructor() {}

    /**
     * An S2Cell always corresponds to a particular S2CellId. The other
     * constructors are just convenience methods.
     */
    constructor(id: S2CellId) {
        init(id)
    }

    // Convenience methods.
    constructor(p: S2Point) {
        init(fromPoint(p))
    }

    constructor(ll: S2LatLng) {
        init(fromLatLng(ll))
    }

    fun id(): S2CellId {
        return cellId
    }

    fun face(): Int {
        return face.toInt()
    }

    fun level(): Byte {
        return level
    }

    fun orientation(): Byte {
        return orientation
    }

    val isLeaf: Boolean
        get() = level == S2CellId.kMaxLevel.toByte()

    fun getVertex(k: Int): S2Point {
        return normalize(getVertexRaw(k))
    }

    /**
     * Return the k-th vertex of the cell (k = 0,1,2,3). Vertices are returned in
     * CCW order. The points returned by GetVertexRaw are not necessarily unit
     * length.
     */
    fun getVertexRaw(k: Int): S2Point {
        // Vertices are returned in the order SW, SE, NE, NW.
        return S2Projections.faceUvToXyz(face.toInt(), uv[0][k shr 1 xor (k and 1)], uv[1][k shr 1])
    }

    fun getEdge(k: Int): S2Point {
        return normalize(getEdgeRaw(k))
    }

    fun getEdgeRaw(k: Int): S2Point {
        return when (k) {
            0 -> S2Projections.getVNorm(face.toInt(), uv[1][0]) // South
            1 -> S2Projections.getUNorm(face.toInt(), uv[0][1]) // East
            2 -> unaryMinus(S2Projections.getVNorm(face.toInt(), uv[1][1])) // North
            else -> unaryMinus(S2Projections.getUNorm(face.toInt(), uv[0][0])) // West
        }
    }

    /**
     * Return the inward-facing normal of the great circle passing through the
     * edge from vertex k to vertex k+1 (mod 4). The normals returned by
     * GetEdgeRaw are not necessarily unit length.
     *
     *
     * If this is not a leaf cell, set children[0..3] to the four children of
     * this cell (in traversal order) and return true. Otherwise returns false.
     * This method is equivalent to the following:
     *
     *
     * for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
     * children[i] = S2Cell(id);
     *
     *
     * except that it is more than two times faster.
     */
    fun subdivide(children: Array<S2Cell>): Boolean {
        // This function is equivalent to just iterating over the child cell ids
        // and calling the S2Cell constructor, but it is about 2.5 times faster.
        if (cellId.isLeaf()) {
            return false
        }

        // Compute the cell midpoint in uv-space.
        val uvMid = centerUV

        // Create four children with the appropriate bounds.
        var id = cellId.childBegin()
        var pos = 0
        while (pos < 4) {
            val child = children[pos]
            child.face = face
            child.level = (level + 1).toByte()
            child.orientation = (orientation.toInt() xor S2.posToOrientation(pos)).toByte()
            child.cellId = id
            val ij = S2.posToIJ(orientation.toInt(), pos)
            for (d in 0..1) {
                // The dimension 0 index (i/u) is in bit 1 of ij.
                val m = 1 - (ij shr 1 - d and 1)
                child.uv[d][m] = uvMid[d]
                child.uv[d][1 - m] = uv[d][1 - m]
            }
            ++pos
            id = id.next()
        }
        return true
    }

    /**
     * Return the direction vector corresponding to the center in (s,t)-space of
     * the given cell. This is the point at which the cell is divided into four
     * subcells; it is not necessarily the centroid of the cell in (u,v)-space or
     * (x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit
     * length.
     */
    val center: S2Point
        get() = normalize(centerRaw)
    val centerRaw: S2Point
        get() = cellId.toPointRaw()// TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).

    /**
     * Return the center of the cell in (u,v) coordinates (see `S2Projections`). Note that the center of the cell is defined as the point
     * at which it is recursively subdivided into four children; in general, it is
     * not at the midpoint of the (u,v) rectangle covered by the cell
     */
    val centerUV: R2Vector
        get() {
            val (_, i, j) = cellId.toFaceIJOrientation()
            val cellSize = 1 shl S2CellId.kMaxLevel - level

            // TODO(dbeaumont): Figure out a better naming of the variables here (and elsewhere).
            val si = (i and -cellSize) * 2 + cellSize - MAX_CELL_SIZE
            val x = S2Projections.stToUV(1.0 / MAX_CELL_SIZE * si)
            val sj = (j and -cellSize) * 2 + cellSize - MAX_CELL_SIZE
            val y = S2Projections.stToUV(1.0 / MAX_CELL_SIZE * sj)
            return R2Vector(x, y)
        }

    /**
     * Return the approximate area of this cell. This method is accurate to within
     * 3% percent for all cell sizes and accurate to within 0.1% for cells at
     * level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to
     * compute.
     */
    fun approxArea(): Double {

        // All cells at the first two levels have the same area.
        if (level < 2) {
            return averageArea(level.toInt())
        }

        // First, compute the approximate area of the cell when projected
        // perpendicular to its normal. The cross product of its diagonals gives
        // the normal, and the length of the normal is twice the projected area.
        val flatArea = 0.5 * crossProd(
                minus(getVertex(2), getVertex(0)), minus(getVertex(3), getVertex(1))).norm()

        // Now, compensate for the curvature of the cell surface by pretending
        // that the cell is shaped like a spherical cap. The ratio of the
        // area of a spherical cap to the area of its projected disc turns out
        // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
        // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
        // Here we set Pi*r*r == flat_area to find the equivalent disc.
        return flatArea * 2 / (1 + Math.sqrt(1 - Math.min(S2.M_1_PI * flatArea, 1.0)))
    }

    /**
     * Return the area of this cell as accurately as possible. This method is more
     * expensive but it is accurate to 6 digits of precision even for leaf cells
     * (whose area is approximately 1e-18).
     */
    fun exactArea(): Double {
        val v0 = getVertex(0)
        val v1 = getVertex(1)
        val v2 = getVertex(2)
        val v3 = getVertex(3)
        return S2.area(v0, v1, v2) + S2.area(v0, v2, v3)
    }

    // //////////////////////////////////////////////////////////////////////
    // S2Region interface (see {@code S2Region} for details):
    public override fun clone(): S2Region {
        val clone = S2Cell()
        clone.face = face
        clone.level = level
        clone.orientation = orientation
        clone.uv = uv.clone()
        return clone
    }

    // Use the cell center in (u,v)-space as the cap axis. This vector is
    // very close to GetCenter() and faster to compute. Neither one of these
    // vectors yields the bounding cap with minimal surface area, but they
    // are both pretty close.
    //
    // It's possible to show that the two vertices that are furthest from
    // the (u,v)-origin never determine the maximum cap size (this is a
    // possible future optimization).
    override val capBound: S2Cap
        get() {
            // Use the cell center in (u,v)-space as the cap axis. This vector is
            // very close to GetCenter() and faster to compute. Neither one of these
            // vectors yields the bounding cap with minimal surface area, but they
            // are both pretty close.
            //
            // It's possible to show that the two vertices that are furthest from
            // the (u,v)-origin never determine the maximum cap size (this is a
            // possible future optimization).
            val u = 0.5 * (uv[0][0] + uv[0][1])
            val v = 0.5 * (uv[1][0] + uv[1][1])
            var cap = fromCenterHeight(normalize(S2Projections.faceUvToXyz(face.toInt(), u, v)), 0.0)
            for (k in 0..3) {
                cap = cap.addPoint(getVertex(k))
            }
            return cap
        }// Except for cells at level 0, the latitude and longitude extremes are

    // attained at the vertices. Furthermore, the latitude range is
    // determined by one pair of diagonally opposite vertices and the
    // longitude range is determined by the other pair.
    //
    // We first determine which corner (i,j) of the cell has the largest
    // absolute latitude. To maximize latitude, we want to find the point in
    // the cell that has the largest absolute z-coordinate and the smallest
    // absolute x- and y-coordinates. To do this we look at each coordinate
    // (u and v), and determine whether we want to minimize or maximize that
    // coordinate based on the axis direction and the cell's (u,v) quadrant.
    // 35.26 degrees
    override val rectBound: S2LatLngRect
        get() {
            if (level > 0) {
                // Except for cells at level 0, the latitude and longitude extremes are
                // attained at the vertices. Furthermore, the latitude range is
                // determined by one pair of diagonally opposite vertices and the
                // longitude range is determined by the other pair.
                //
                // We first determine which corner (i,j) of the cell has the largest
                // absolute latitude. To maximize latitude, we want to find the point in
                // the cell that has the largest absolute z-coordinate and the smallest
                // absolute x- and y-coordinates. To do this we look at each coordinate
                // (u and v), and determine whether we want to minimize or maximize that
                // coordinate based on the axis direction and the cell's (u,v) quadrant.
                val u = uv[0][0] + uv[0][1]
                val v = uv[1][0] + uv[1][1]
                val i = if (S2Projections.getUAxis(face.toInt()).z() == 0.0) if (u < 0) 1 else 0 else if (u > 0) 1 else 0
                val j = if (S2Projections.getVAxis(face.toInt()).z() == 0.0) if (v < 0) 1 else 0 else if (v > 0) 1 else 0
                var lat = R1Interval.fromPointPair(getLatitude(i, j), getLatitude(1 - i, 1 - j))
                lat = lat.expanded(MAX_ERROR).intersection(fullLat())
                if (lat.lo == -S2.M_PI_2 || lat.hi == S2.M_PI_2) {
                    return S2LatLngRect(lat, full)
                }
                val lng = S1Interval.fromPointPair(getLongitude(i, 1 - j), getLongitude(1 - i, j))
                return S2LatLngRect(lat, lng.expanded(MAX_ERROR))
            }
            return when (face.toInt()) {
                0 -> S2LatLngRect(
                        R1Interval(-S2.M_PI_4, S2.M_PI_4), S1Interval(-S2.M_PI_4, S2.M_PI_4))
                1 -> S2LatLngRect(
                        R1Interval(-S2.M_PI_4, S2.M_PI_4), S1Interval(S2.M_PI_4, 3 * S2.M_PI_4))
                2 -> S2LatLngRect(
                        R1Interval(POLE_MIN_LAT, S2.M_PI_2), S1Interval(-S2.M_PI, S2.M_PI))
                3 -> S2LatLngRect(
                        R1Interval(-S2.M_PI_4, S2.M_PI_4), S1Interval(3 * S2.M_PI_4, -3 * S2.M_PI_4))
                4 -> S2LatLngRect(
                        R1Interval(-S2.M_PI_4, S2.M_PI_4), S1Interval(-3 * S2.M_PI_4, -S2.M_PI_4))
                else -> S2LatLngRect(
                        R1Interval(-S2.M_PI_2, -POLE_MIN_LAT), S1Interval(-S2.M_PI, S2.M_PI))
            }
        }

    override fun mayIntersect(cell: S2Cell): Boolean {
        return cellId.intersects(cell.cellId)
    }

    override operator fun contains(p: S2Point): Boolean {
        // We can't just call XYZtoFaceUV, because for points that lie on the
        // boundary between two faces (i.e. u or v is +1/-1) we need to return
        // true for both adjacent cells.
        val uvPoint = S2Projections.faceXyzToUv(face.toInt(), p) ?: return false
        return uvPoint.x() >= uv[0][0] && uvPoint.x() <= uv[0][1] && uvPoint.y() >= uv[1][0] && uvPoint.y() <= uv[1][1]
    }

    // The point 'p' does not need to be normalized.
    override fun contains(cell: S2Cell): Boolean {
        return cellId.contains(cell.cellId)
    }

    private fun init(id: S2CellId) {
        cellId = id
        val (face1, i, j, orientation1) = id.toFaceIJOrientation(true)
        face = face1.toByte()
        orientation = orientation1!!.toByte() // Compress int to a byte.
        level = id.level().toByte()
        val cellSize = 1 shl S2CellId.kMaxLevel - level
        val ij = intArrayOf(i, j)
        for (d in 0..1) {
            // Compute the cell bounds in scaled (i,j) coordinates.
            val sijLo = (ij[d] and -cellSize) * 2 - MAX_CELL_SIZE
            val sijHi = sijLo + cellSize * 2
            uv[d][0] = S2Projections.stToUV(1.0 / MAX_CELL_SIZE * sijLo)
            uv[d][1] = S2Projections.stToUV(1.0 / MAX_CELL_SIZE * sijHi)
        }
    }

    // Internal method that does the actual work in the constructors.
    private fun getLatitude(i: Int, j: Int): Double {
        val p = S2Projections.faceUvToXyz(face.toInt(), uv[0][i], uv[1][j])
        return Math.atan2(p.z(), Math.sqrt(p.x() * p.x() + p.y() * p.y()))
    }

    private fun getLongitude(i: Int, j: Int): Double {
        val p = S2Projections.faceUvToXyz(face.toInt(), uv[0][i], uv[1][j])
        return Math.atan2(p.y(), p.x())
    }

    // Return the latitude or longitude of the cell vertex given by (i,j),
    // where "i" and "j" are either 0 or 1.
    override fun toString(): String {
        return "[$face, $level, $orientation, $cellId]"
    }

    override fun hashCode(): Int {
        var value = 17
        value = 37 * (37 * (37 * value + face) + orientation) + level
        return 37 * value + id().hashCode()
    }

    override fun equals(that: Any?): Boolean {
        if (that is S2Cell) {
            val thatCell = that
            return face == thatCell.face && level == thatCell.level && orientation == thatCell.orientation && cellId!!.equals(thatCell.cellId)
        }
        return false
    }

    fun averageArea(): Double = averageArea(level.toInt())

    companion object {
        private const val MAX_CELL_SIZE = 1 shl S2CellId.kMaxLevel

        // This is a static method in order to provide named parameters.
        @JvmStatic
        fun fromFacePosLevel(face: Int, pos: Byte, level: Int): S2Cell {
            return S2Cell(S2CellId.fromFacePosLevel(face, pos.toULong(), level))
        }

        fun fromFace(face: Int): S2Cell {
            return S2Cell(S2CellId.fromFace(face))
        }
        /**
         * Return the average area for cells at the given level.
         */
        /**
         * Return the average area of cells at this level. This is accurate to within
         * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
         * compute.
         */
        @JvmStatic
        fun averageArea(level: Int): Double {
            return S2Projections.AVG_AREA.getValue(level)
        }

        // We grow the bounds slightly to make sure that the bounding rectangle
        // also contains the normalized versions of the vertices. Note that the
        // maximum result magnitude is Pi, with a floating-point exponent of 1.
        // Therefore adding or subtracting 2**-51 will always change the result.
        private const val MAX_ERROR = 1.0 / (1L shl 51)

        // The 4 cells around the equator extend to +/-45 degrees latitude at the
        // midpoints of their top and bottom edges. The two cells covering the
        // poles extend down to +/-35.26 degrees at their vertices.
        // adding kMaxError (as opposed to the C version) because of asin and atan2
        // roundoff errors
        private val POLE_MIN_LAT = Math.asin(Math.sqrt(1.0 / 3.0)) - MAX_ERROR
    }
}