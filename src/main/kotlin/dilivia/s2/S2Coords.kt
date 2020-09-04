package dilivia.s2

import com.google.common.geometry.S2.*
import dilivia.s2.math.R2Point
import kotlin.math.*

// S2 is a namespace for constants and simple utility functions that are used
// throughout the S2 library.  The name "S2" is derived from the mathematical
// symbol for the two-dimensional unit sphere (note that the "2" refers to the
// dimension of the surface, not the space it is embedded in).


// We have implemented three different projections from cell-space (s,t) to
// cube-space (u,v): linear, quadratic, and tangent.  They have the following
// tradeoffs:
//
//   Linear - This is the fastest transformation, but also produces the least
//   uniform cell sizes.  Cell areas vary by a factor of about 5.2, with the
//   largest cells at the center of each face and the smallest cells in
//   the corners.
//
//   Tangent - Transforming the coordinates via atan() makes the cell sizes
//   more uniform.  The areas vary by a maximum ratio of 1.4 as opposed to a
//   maximum ratio of 5.2.  However, each call to atan() is about as expensive
//   as all of the other calculations combined when converting from points to
//   cell ids, i.e. it reduces performance by a factor of 3.
//
//   Quadratic - This is an approximation of the tangent projection that
//   is much faster and produces cells that are almost as uniform in size.
//   It is about 3 times faster than the tangent projection for converting
//   cell ids to points or vice versa.  Cell areas vary by a maximum ratio of
//   about 2.1.
//
// Here is a table comparing the cell uniformity using each projection.  "Area
// ratio" is the maximum ratio over all subdivision levels of the largest cell
// area to the smallest cell area at that level, "edge ratio" is the maximum
// ratio of the longest edge of any cell to the shortest edge of any cell at
// the same level, and "diag ratio" is the ratio of the longest diagonal of
// any cell to the shortest diagonal of any cell at the same level.  "ToPoint"
// and "FromPoint" are the times in microseconds required to convert cell ids
// to and from points (unit vectors) respectively.  "ToPointRaw" is the time
// to convert to a non-unit-length vector, which is all that is needed for
// some purposes.
//
//               Area    Edge    Diag   ToPointRaw  ToPoint  FromPoint
//              Ratio   Ratio   Ratio             (microseconds)
// -------------------------------------------------------------------
// Linear:      5.200   2.117   2.959      0.020     0.087     0.085
// Tangent:     1.414   1.414   1.704      0.237     0.299     0.258
// Quadratic:   2.082   1.802   1.932      0.033     0.096     0.108
//
// The worst-case cell aspect ratios are about the same with all three
// projections.  The maximum ratio of the longest edge to the shortest edge
// within the same cell is about 1.4 and the maximum ratio of the diagonals
// within the same cell is about 1.7.
//
// This data was produced using s2cell_test and s2cell_id_test.

enum class Projections {
    LINEAR, TAN, QUADRATIC
}

interface S2Projection {
    // Convert an s- or t-value to the corresponding u- or v-value.  This is
    // a non-linear transformation from [-1,1] to [-1,1] that attempts to
    // make the cell sizes more uniform.
    fun stToUv(s: Double): Double

    // The inverse of the STtoUV transformation.  Note that it is not always
    // true that UVtoST(STtoUV(x)) == x due to numerical errors.
    fun uvToSt(u: Double): Double
}

object S2LinearProjection : S2Projection {
    override fun stToUv(s: Double): Double = 2 * s - 1
    override fun uvToSt(u: Double): Double = 0.5 * (u + 1)
}

object S2TanProjection : S2Projection {

    override fun stToUv(s: Double): Double {
        // Unfortunately, tan(M_PI_4) is slightly less than 1.0.  This isn't due to
        // a flaw in the implementation of tan(), it's because the derivative of
        // tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
        // point numbers on either side of the infinite-precision value of pi/4 have
        // tangents that are slightly below and slightly above 1.0 when rounded to
        // the nearest double-precision result.

        val t = tan(M_PI_2 * s - M_PI_4)
        return t + (1.0 / (1L shl 53)) * t
    }

    override fun uvToSt(u: Double): Double {
        val a = atan(u)
        return (2 * M_1_PI) * (a + M_PI_4)
    }

}

object S2QuadraticProjection : S2Projection {
    override fun stToUv(s: Double): Double = if (s >= 0.5) (1.0 / 3.0) * (4 * s * s - 1)
    else (1.0 / 3.0) * (1 - 4 * (1 - s) * (1 - s))

    override fun uvToSt(u: Double): Double = if (u >= 0) 0.5 * sqrt(1 + 3 * u)
    else 1 - 0.5 * sqrt(1 - 3 * u)
}

@ExperimentalUnsignedTypes
object S2Coords {

    // The canonical Hilbert traversal order looks like an inverted 'U':
    // the subcells are visited in the order (0,0), (0,1), (1,1), (1,0).
    // The following tables encode the traversal order for various
    // orientations of the Hilbert curve (axes swapped and/or directions
    // of the axes reversed).

    // Together these flags define a cell orientation.  If 'kSwapMask'
    // is true, then canonical traversal order is flipped around the
    // diagonal (i.e. i and j are swapped with each other).  If
    // 'kInvertMask' is true, then the traversal order is rotated by 180
    // degrees (i.e. the bits of i and j are inverted, or equivalently,
    // the axis directions are reversed).
    internal val kSwapMask = 0x01
    internal val kInvertMask = 0x02

    // kIJtoPos[orientation][ij] -> pos
    //
    // Given a cell orientation and the (i,j)-index of a subcell (0=(0,0),
    // 1=(0,1), 2=(1,0), 3=(1,1)), return the order in which this subcell is
    // visited by the Hilbert curve (a position in the range [0..3]).
    internal val kIJtoPos = arrayOf(
            //        (0,0) (0,1) (1,0) (1,1)
            intArrayOf(0,    1,    3,    2),  // canonical order
            intArrayOf(0,    3,    1,    2),  // axes swapped
            intArrayOf(2,    3,    1,    0),  // bits inverted
            intArrayOf(2,    1,    3,    0),  // swapped & inverted
    ) //[4][4];

    // kPosToIJ[orientation][pos] -> ij
    //
    // Return the (i,j) index of the subcell at the given position 'pos' in the
    // Hilbert curve traversal order with the given orientation.  This is the
    // inverse of the previous table:
    //
    //   kPosToIJ[r][kIJtoPos[r][ij]] == ij
    internal val kPosToIJ = arrayOf(
            //         0  1  2  3
            intArrayOf(0, 1, 3, 2),    // canonical order:    (0,0), (0,1), (1,1), (1,0)
            intArrayOf(0, 2, 3, 1),    // axes swapped:       (0,0), (1,0), (1,1), (0,1)
            intArrayOf(3, 2, 0, 1),    // bits inverted:      (1,1), (1,0), (0,0), (0,1)
            intArrayOf(3, 1, 0, 2),    // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
    ) //[4][4];

    // kPosToOrientation[pos] -> orientation_modifier
    //
    // Return a modifier indicating how the orientation of the child subcell
    // with the given traversal position [0..3] is related to the orientation
    // of the parent cell.  The modifier should be XOR-ed with the parent
    // orientation to obtain the curve orientation in the child.
    private val kPosToOrientation = intArrayOf(
        kSwapMask,
        0,
        0,
        kInvertMask + kSwapMask,
    )

    // The U,V,W axes for each face.
    val kFaceUVWAxes //[6][3][3]
    = arrayOf(
        arrayOf(
            intArrayOf( 0,  1,  0 ),
            intArrayOf( 0,  0,  1 ),
            intArrayOf( 1,  0,  0 )
        ),
        arrayOf(
            intArrayOf(-1,  0,  0 ),
            intArrayOf( 0,  0,  1 ),
            intArrayOf( 0,  1,  0 )
        ),
        arrayOf(
            intArrayOf(-1,  0,  0 ),
            intArrayOf( 0, -1,  0 ),
            intArrayOf( 0,  0,  1 )
        ),
        arrayOf(
            intArrayOf( 0,  0, -1 ),
            intArrayOf( 0, -1,  0 ),
            intArrayOf(-1,  0,  0 )
        ),
        arrayOf(
            intArrayOf( 0,  0, -1 ),
            intArrayOf( 1,  0,  0 ),
            intArrayOf( 0, -1,  0 )
        ),
        arrayOf(
            intArrayOf( 0,  1,  0 ),
            intArrayOf( 1,  0,  0 ),
            intArrayOf( 0,  0, -1 )
        )
    )

    // The precomputed neighbors of each face (see GetUVWFace).
    val kFaceUVWFaces //[6][3][2]
    = arrayOf(
        arrayOf( intArrayOf( 4, 1 ), intArrayOf( 5, 2 ), intArrayOf( 3, 0 ) ),
        arrayOf( intArrayOf( 0, 3 ), intArrayOf( 5, 2 ), intArrayOf( 4, 1 ) ),
        arrayOf( intArrayOf( 0, 3 ), intArrayOf( 1, 4 ), intArrayOf( 5, 2 ) ),
        arrayOf( intArrayOf( 2, 5 ), intArrayOf( 1, 4 ), intArrayOf( 0, 3 ) ),
        arrayOf( intArrayOf( 2, 5 ), intArrayOf( 3, 0 ), intArrayOf( 1, 4 ) ),
        arrayOf( intArrayOf( 4, 1 ), intArrayOf( 3, 0 ), intArrayOf( 2, 5 ) )
    )

    @JvmField
    val S2_PROJECTION = Projections.QUADRATIC

    var projection: S2Projection = S2QuadraticProjection

    // This is the number of levels needed to specify a leaf cell.  This
    // constant is defined here so that the S2::Metric class and the conversion
    // functions below can be implemented without including s2cell_id.h.  Please
    // see s2cell_id.h for other useful constants and conversion functions.
    const val kMaxCellLevel = 30

    // The maximum index of a valid leaf cell plus one.  The range of valid leaf
    // cell indices is [0..kLimitIJ-1].
    val kLimitIJ = 1 shl kMaxCellLevel  // == S2CellId::kMaxSize

    // The maximum value of an si- or ti-coordinate.  The range of valid (si,ti)
    // values is [0..kMaxSiTi].
    val kMaxSiTi = 1.toUInt() shl (kMaxCellLevel + 1)

    // Convert an s- or t-value to the corresponding u- or v-value.  This is
    // a non-linear transformation from [-1,1] to [-1,1] that attempts to
    // make the cell sizes more uniform.
    fun stToUv(s: Double): Double = projection.stToUv(s)

    // The inverse of the STtoUV transformation.  Note that it is not always
    // true that UVtoST(STtoUV(x)) == x due to numerical errors.
    fun uvToSt(u: Double): Double = projection.uvToSt(u)

    // Convert the i- or j-index of a leaf cell to the minimum corresponding s-
    // or t-value contained by that cell.  The argument must be in the range
    // [0..2**30], i.e. up to one position beyond the normal range of valid leaf
    // cell indices.
    fun ijToStMin(i: Int): Double {
        assert(i in 0..kLimitIJ)
        return (1.0 / kLimitIJ) * i
    }

    // Return the i- or j-index of the leaf cell containing the given
    // s- or t-value.  If the argument is outside the range spanned by valid
    // leaf cell indices, return the index of the closest valid leaf cell (i.e.,
    // return values are clamped to the range of valid leaf cell indices).
    fun stToIj(s: Double): Int = max(0, min(kLimitIJ - 1, (kLimitIJ * s - 0.5).roundToInt()))

    // Convert an si- or ti-value to the corresponding s- or t-value.
    fun siTiToSt(si: UInt): Double {
        assert(si <= kMaxSiTi)
        return (1.0 / kMaxSiTi.toDouble()) * si.toDouble()
    }

    // Return the si- or ti-coordinate that is nearest to the given s- or
    // t-value.  The result may be outside the range of valid (si,ti)-values.
    fun stToSiTi(s: Double): UInt {
        // kMaxSiTi == 2^31, so the result doesn't fit in an int32 when s == 1.
        return (s * kMaxSiTi.toDouble()).roundToLong().toUInt()
    }

    // Convert (face, u, v) coordinates to a direction vector (not
    // necessarily unit length).
    fun faceUVtoXYZ(faceUV: FaceUV): S2Point = faceUVtoXYZ(faceUV.face, faceUV.u, faceUV.v)
    fun faceUVtoXYZ(face: Int, u: Double, v: Double): S2Point {
        return when (face) {
            0 -> S2Point(1.0, u, v)
            1 -> S2Point(-u, 1.0, v)
            2 -> S2Point(-u, -v, 1.0)
            3 -> S2Point(-1.0, -v, -u)
            4 -> S2Point(v, -1.0, -u)
            else -> S2Point(v, u, -1.0)
        }
    }

    // If the dot product of p with the given face normal is positive,
    // set the corresponding u and v values (which may lie outside the range
    // [-1,1]) and return true.  Otherwise return false.
    fun faceXYZtoUV(face: Int, p: S2Point): R2Point? {
        if (face < 3) {
            if (p[face] <= 0) return null
        } else {
            if (p[face - 3] >= 0) return null
        }
        return validFaceXYZtoUV(face, p)
    }

    // Given a *valid* face for the given point p (meaning that dot product
    // of p with the face normal is positive), return the corresponding
    // u and v values (which may lie outside the range [-1,1]).
    fun validFaceXYZtoUV(face: Int, p: S2Point): R2Point {
        Assertions.assert { p.dotProd(getNorm(face)) > 0 }
        return when (face) {
            0 -> R2Point(x = p[1] / p[0], y = p[2] / p[0])
            1 -> R2Point(x = -p[0] / p[1], y = p[2] / p[1])
            2 -> R2Point(x = -p[0] / p[2], y = -p[1] / p[2])
            3 -> R2Point(x = p[2] / p[0], y = p[1] / p[0])
            4 -> R2Point(x = p[2] / p[1], y = -p[0] / p[1])
            else -> R2Point(x = -p[1] / p[2], y = -p[0] / p[2])
        }
    }

    // Transform the given point P to the (u,v,w) coordinate frame of the given
    // face (where the w-axis represents the face normal).
    fun faceXYZtoUVW(face: Int, p: S2Point): S2Point {
        // The result coordinates are simply the dot products of P with the (u,v,w)
        // axes for the given face (see kFaceUVWAxes).
        return when (face) {
            0 -> S2Point( p.y(),  p.z(),  p.x())
            1 -> S2Point(-p.x(),  p.z(),  p.y())
            2 -> S2Point(-p.x(), -p.y(),  p.z())
            3 -> S2Point(-p.z(), -p.y(), -p.x())
            4 -> S2Point(-p.z(),  p.x(), -p.y())
            else -> S2Point( p.y(),  p.x(), -p.z())
        }
    }

    // Return the face containing the given direction vector.  (For points on
    // the boundary between faces, the result is arbitrary but repeatable.)
    fun getFace(p: S2Point): Int {
        var face = p.largestAbsComponent()
        if (p[face] < 0) face += 3
        return face;
    }

    // Convert a direction vector (not necessarily unit length) to
    // (face, u, v) coordinates.
    fun xyzToFaceUV(p: S2Point): FaceUV {
        val face = getFace(p)
        val uv = validFaceXYZtoUV(face, p)
        return FaceUV(face = face, u = uv.x(), v = uv.y())
    }

    // Convert a direction vector (not necessarily unit length) to
    // (face, si, ti) coordinates and, if p is exactly equal to the center of a
    // cell, return the level of this cell (-1 otherwise).
    fun xyzToFaceSiTi(p: S2Point): Pair<Int, FaceSiTi> {
        val (face, u, v) = xyzToFaceUV(p)
        val si = stToSiTi(projection.uvToSt(u))
        val ti = stToSiTi(projection.uvToSt(v))
        val faceSiTi = FaceSiTi(face, si, ti)
        // If the levels corresponding to si,ti are not equal, then p is not a cell
        // center.  The si,ti values 0 and kMaxSiTi need to be handled specially
        // because they do not correspond to cell centers at any valid level; they
        // are mapped to level -1 by the code below.
        val level = kMaxCellLevel - (si or kMaxSiTi).countTrailingZeroBits()
        if (level < 0 || level != kMaxCellLevel - (ti or kMaxSiTi).countTrailingZeroBits()) {
            return -1 to faceSiTi
        }
        assert(level <= kMaxCellLevel);
        // In infinite precision, this test could be changed to ST == SiTi. However,
        // due to rounding errors, UVtoST(XYZtoFaceUV(FaceUVtoXYZ(STtoUV(...)))) is
        // not idempotent. On the other hand, center_raw is computed exactly the same
        // way p was originally computed (if it is indeed the center of an S2Cell):
        // the comparison can be exact.
        val center = faceSiTitoXYZ(faceSiTi).normalize()
        return (if(p == center) level else -1) to faceSiTi
    }

    // Convert (face, si, ti) coordinates to a direction vector (not necessarily
    // unit length).
    fun faceSiTitoXYZ(faceSiTi: FaceSiTi): S2Point = faceSiTitoXYZ(faceSiTi.face, faceSiTi.si, faceSiTi.ti)
    fun faceSiTitoXYZ(face: Int, si: UInt, ti: UInt): S2Point {
        val u = projection.stToUv(siTiToSt(si))
        val v = projection.stToUv(siTiToSt(ti))
        return faceUVtoXYZ(face, u, v)
    }

    // Return the right-handed normal (not necessarily unit length) for an
    // edge in the direction of the positive v-axis at the given u-value on
    // the given face.  (This vector is perpendicular to the plane through
    // the sphere origin that contains the given edge.)
    fun getUNorm(face: Int, u: Double): S2Point = when (face) {
        0 -> S2Point(u, -1.0, 0.0)
        1 -> S2Point(1.0, u, 0.0)
        2 -> S2Point(1.0, 0.0, u)
        3 -> S2Point(-u, 0.0, 1.0)
        4 -> S2Point(0.0, -u, 1.0)
        else -> S2Point(0.0, -1.0, -u)
    }

    // Return the right-handed normal (not necessarily unit length) for an
    // edge in the direction of the positive u-axis at the given v-value on
    // the given face.
    fun getVNorm(face: Int, v: Double): S2Point {
        return when (face) {
            0 -> S2Point(-v, 0.0, 1.0)
            1 -> S2Point(0.0, -v, 1.0)
            2 -> S2Point(0.0, -1.0, -v)
            3 -> S2Point(v, -1.0, 0.0)
            4 -> S2Point(1.0, v, 0.0)
            else -> S2Point(1.0, 0.0, v)
        }
    }

    // Return the unit-length normal, u-axis, or v-axis for the given face.
    fun getNorm(face: Int): S2Point = getUVWAxis(face, 2)
    fun getUAxis(face: Int): S2Point = getUVWAxis(face, 0)
    fun getVAxis(face: Int): S2Point = getUVWAxis(face, 1)

    // Return the given axis of the given face (u=0, v=1, w=2).
    fun getUVWAxis(face: Int, axis: Int): S2Point {
        val p = kFaceUVWAxes[face][axis]
        return S2Point(p[0], p[1], p[2])
    }

    // With respect to the (u,v,w) coordinate system of a given face, return the
    // face that lies in the given direction (negative=0, positive=1) of the
    // given axis (u=0, v=1, w=2).  For example, GetUVWFace(4, 0, 1) returns the
    // face that is adjacent to face 4 in the positive u-axis direction.
    fun getUVWFace(face: Int, axis: Int, direction: Int): Int {
        assert(face in 0..5)
        assert(axis in 0..2)
        assert(direction in 0..1)
        return kFaceUVWFaces[face][axis][direction]
    }

}

data class FaceUV(val face: Int, val u: Double, val v: Double)
data class FaceSiTi(val face: Int, val si: UInt, val ti: UInt)