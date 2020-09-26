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
    fun stToUV(s: Double): Double

    // The inverse of the STtoUV transformation.  Note that it is not always
    // true that UVtoST(STtoUV(x)) == x due to numerical errors.
    fun uvToST(u: Double): Double

    val kMinAngleSpan: LengthMetric
    val kMaxAngleSpan: LengthMetric
    val kAvgAngleSpan: LengthMetric

    val kMinWidth: LengthMetric
    val kMaxWidth: LengthMetric
    val kAvgWidth: LengthMetric

    val kMinEdge: LengthMetric
    val kMaxEdge: LengthMetric
    val kAvgEdge: LengthMetric
    
    val kMinDiag: LengthMetric
    val kMaxDiag: LengthMetric
    val kAvgDiag: LengthMetric
    
    val kMinArea: AreaMetric
    val kMaxArea: AreaMetric
    val kAvgArea: AreaMetric

    val kMaxEdgeAspect: Double
    val kMaxDiagAspect: Double

}

object S2LinearProjection : S2Projection {
    override fun stToUV(s: Double): Double = 2 * s - 1
    override fun uvToST(u: Double): Double = 0.5 * (u + 1)
    
    override val kMinAngleSpan: LengthMetric = LengthMetric(1.0)                                          // 1.000 
    override val kMaxAngleSpan: LengthMetric = LengthMetric(2.0)                                          // 2.000
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                                     // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(sqrt(2.0 / 3.0))                                     // 0.816
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.411459345844456965)                             // 1.411

    override val kMinEdge: LengthMetric = LengthMetric(2 * sqrt(2.0) / 3.0)                            // 0.943
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.440034192955603643)                              // 1.440

    override val kMinDiag: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)                          // 0.943
    override val kMaxDiag: LengthMetric = LengthMetric(2.0 * sqrt(2.0))                                // 2.828
    override val kAvgDiag: LengthMetric = LengthMetric(2.031817866418812674)                              // 2.032
    
    override val kMinArea: AreaMetric = AreaMetric(4.0 / (3.0 * sqrt(3.0)))                            // 0.770
    override val kMaxArea: AreaMetric = AreaMetric(4.0)                                                   // 4.000
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6);                                         // 2.094

    override val kMaxEdgeAspect: Double = sqrt(2.0)                                                          // 1.414
    override val kMaxDiagAspect: Double = sqrt(3.0)                                                          // 1.732

}

object S2TanProjection : S2Projection {

    override fun stToUV(s: Double): Double {
        // Unfortunately, tan(M_PI_4) is slightly less than 1.0.  This isn't due to
        // a flaw in the implementation of tan(), it's because the derivative of
        // tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
        // point numbers on either side of the infinite-precision value of pi/4 have
        // tangents that are slightly below and slightly above 1.0 when rounded to
        // the nearest double-precision result.

        val t = tan(M_PI_2 * s - M_PI_4)
        return t + (1.0 / (1L shl 53)) * t
    }

    override fun uvToST(u: Double): Double {
        val a = atan(u)
        return (2 * M_1_PI) * (a + M_PI_4)
    }
    
    override val kMinAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                                     // 1.571
    override val kMaxAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                                     // 1.571
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                                     // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(M_PI / (2 * sqrt(2.0)))                        // 1.111
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.437318638925160885)                             // 1.437

    override val kMinEdge: LengthMetric = LengthMetric(M_PI / (2 * sqrt(2.0)))                         // 1.111
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.461667032546739266)                              // 1.462

    override val kMinDiag: LengthMetric = LengthMetric(M_PI * sqrt(2.0) / 3)                           // 1.481
    override val kMaxDiag: LengthMetric = LengthMetric(M_PI * sqrt(2.0 / 3.0))                         // 2.565
    override val kAvgDiag: LengthMetric = LengthMetric(2.063623197195635753)                              // 2.064

    override val kMinArea: AreaMetric = AreaMetric((M_PI*M_PI) / (4.0 * sqrt(2.0)))                    // 1.745
    override val kMaxArea: AreaMetric = AreaMetric(M_PI * M_PI / 4.0)                                     // 2.467
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6);                                         // 2.094

    override val kMaxEdgeAspect: Double = sqrt(2.0)                                                          // 1.414
    override val kMaxDiagAspect: Double = sqrt(3.0)                                                          // 1.732

}

object S2QuadraticProjection : S2Projection {
    override fun stToUV(s: Double): Double = if (s >= 0.5) (1.0 / 3.0) * (4 * s * s - 1)
    else (1.0 / 3.0) * (1 - 4 * (1 - s) * (1 - s))

    override fun uvToST(u: Double): Double = if (u >= 0) 0.5 * sqrt(1 + 3 * u)
    else 1 - 0.5 * sqrt(1 - 3 * u)
    
    override val kMinAngleSpan: LengthMetric = LengthMetric(4.0 / 3.0)                                    // 1.333
    override val kMaxAngleSpan: LengthMetric = LengthMetric(1.704897179199218452)                         // 1.705
    override val kAvgAngleSpan: LengthMetric = LengthMetric(M_PI / 2)                                     // 1.571

    override val kMinWidth: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)                         // 0.943
    override val kMaxWidth: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgWidth: LengthMetric = LengthMetric(1.434523672886099389)                             // 1.435

    override val kMinEdge: LengthMetric = LengthMetric(2.0 * sqrt(2.0) / 3.0)                          // 0.943
    override val kMaxEdge: LengthMetric = LengthMetric(kMaxAngleSpan.deriv)
    override val kAvgEdge: LengthMetric = LengthMetric(1.459213746386106062)                              // 1.459

    override val kMinDiag: LengthMetric = LengthMetric(8.0 * sqrt(2.0) / 9.0)                          // 1.257
    override val kMaxDiag: LengthMetric = LengthMetric(2.438654594434021032)                              // 2.439
    override val kAvgDiag: LengthMetric = LengthMetric(2.060422738998471683)                              // 2.060

    override val kMinArea: AreaMetric = AreaMetric(8.0 * sqrt(2.0) / 9.0)                              // 1.257
    override val kMaxArea: AreaMetric = AreaMetric(2.635799256963161491)                                  // 2.636
    override val kAvgArea: AreaMetric = AreaMetric(4 * M_PI / 6)                                          // 2.094

    override val kMaxEdgeAspect: Double = 1.442615274452682920                                                  // 1.443
    override val kMaxDiagAspect: Double = sqrt(3.0)                                                          // 1.732

}

//
// This file contains documentation of the various coordinate systems used
// throughout the library.  Most importantly, S2 defines a framework for
// decomposing the unit sphere into a hierarchy of "cells".  Each cell is a
// quadrilateral bounded by four geodesics.  The top level of the hierarchy is
// obtained by projecting the six faces of a cube onto the unit sphere, and
// lower levels are obtained by subdividing each cell into four children
// recursively.  Cells are numbered such that sequentially increasing cells
// follow a continuous space-filling curve over the entire sphere.  The
// transformation is designed to make the cells at each level fairly uniform
// in size.
//
//
////////////////////////// S2Cell Decomposition /////////////////////////
//
// The following methods define the cube-to-sphere projection used by
// the S2Cell decomposition.
//
// In the process of converting a latitude-longitude pair to a 64-bit cell
// id, the following coordinate systems are used:
//
//  (id)
//    An S2CellId is a 64-bit encoding of a face and a Hilbert curve position
//    on that face.  The Hilbert curve position implicitly encodes both the
//    position of a cell and its subdivision level (see s2cell_id.h).
//
//  (face, i, j)
//    Leaf-cell coordinates.  "i" and "j" are integers in the range
//    [0,(2**30)-1] that identify a particular leaf cell on the given face.
//    The (i, j) coordinate system is right-handed on each face, and the
//    faces are oriented such that Hilbert curves connect continuously from
//    one face to the next.
//
//  (face, s, t)
//    Cell-space coordinates.  "s" and "t" are real numbers in the range
//    [0,1] that identify a point on the given face.  For example, the point
//    (s, t) = (0.5, 0.5) corresponds to the center of the top-level face
//    cell.  This point is also a vertex of exactly four cells at each
//    subdivision level greater than zero.
//
//  (face, si, ti)
//    Discrete cell-space coordinates.  These are obtained by multiplying
//    "s" and "t" by 2**31 and rounding to the nearest unsigned integer.
//    Discrete coordinates lie in the range [0,2**31].  This coordinate
//    system can represent the edge and center positions of all cells with
//    no loss of precision (including non-leaf cells).  In binary, each
//    coordinate of a level-k cell center ends with a 1 followed by
//    (30 - k) 0s.  The coordinates of its edges end with (at least)
//    (31 - k) 0s.
//
//  (face, u, v)
//    Cube-space coordinates in the range [-1,1].  To make the cells at each
//    level more uniform in size after they are projected onto the sphere,
//    we apply a nonlinear transformation of the form u=f(s), v=f(t).
//    The (u, v) coordinates after this transformation give the actual
//    coordinates on the cube face (modulo some 90 degree rotations) before
//    it is projected onto the unit sphere.
//
//  (face, u, v, w)
//    Per-face coordinate frame.  This is an extension of the (face, u, v)
//    cube-space coordinates that adds a third axis "w" in the direction of
//    the face normal.  It is always a right-handed 3D coordinate system.
//    Cube-space coordinates can be converted to this frame by setting w=1,
//    while (u,v,w) coordinates can be projected onto the cube face by
//    dividing by w, i.e. (face, u/w, v/w).
//
//  (x, y, z)
//    Direction vector (S2Point).  Direction vectors are not necessarily unit
//    length, and are often chosen to be points on the biunit cube
//    [-1,+1]x[-1,+1]x[-1,+1].  They can be be normalized to obtain the
//    corresponding point on the unit sphere.
//
//  (lat, lng)
//    Latitude and longitude (S2LatLng).  Latitudes must be between -90 and
//    90 degrees inclusive, and longitudes must be between -180 and 180
//    degrees inclusive.
//
// Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are
// right-handed on all six faces.
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
    val kPosToOrientation = intArrayOf(
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
    // functions below can be implemented without including done.s2cell_id.h.  Please
    // see done.s2cell_id.h for other useful constants and conversion functions.
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
    fun stToUv(s: Double): Double = projection.stToUV(s)

    // The inverse of the STtoUV transformation.  Note that it is not always
    // true that UVtoST(STtoUV(x)) == x due to numerical errors.
    fun uvToSt(u: Double): Double = projection.uvToST(u)

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
    fun stToIJ(s: Double): Int = max(0, min(kLimitIJ - 1, (kLimitIJ * s - 0.5).roundToInt()))

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
    fun faceUVtoXYZ(face: Int, uv: R2Point): S2Point = faceUVtoXYZ(face, uv[0], uv[1])
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
    fun validFaceXYZtoUV(face: Int, p: S2Point, uv: MutableR2Point) {
        Assertions.assert { p.dotProd(getNorm(face)) > 0 }
        when (face) {
            0 -> { uv[0] = p[1] / p[0]; uv[1] = p[2] / p[0]; }
            1 -> { uv[0] = -p[0] / p[1]; uv[1] = p[2] / p[1]; }
            2 -> { uv[0] = -p[0] / p[2]; uv[1] = -p[1] / p[2]; }
            3 -> { uv[0] = p[2] / p[0]; uv[1] = p[1] / p[0]; }
            4 -> { uv[0] = p[2] / p[1]; uv[1] = -p[0] / p[1]; }
            else -> { uv[0] = -p[1] / p[2]; uv[1] = -p[0] / p[2]; }
        }
    }
    fun validFaceXYZtoUV(face: Int, p: S2Point): R2Point {
        val uv = MutableR2Point()
        validFaceXYZtoUV(face, p, uv)
        return uv
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
        val si = stToSiTi(projection.uvToST(u))
        val ti = stToSiTi(projection.uvToST(v))
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
        val u = projection.stToUV(siTiToSt(si))
        val v = projection.stToUV(siTiToSt(ti))
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

    fun stToUV(s: Double): Double = projection.stToUV(s)
    fun uvToST(u: Double): Double = projection.uvToST(u)

}

