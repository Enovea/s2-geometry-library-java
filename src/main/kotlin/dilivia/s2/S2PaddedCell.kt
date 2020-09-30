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

import com.google.common.geometry.S2.DBL_EPSILON
import dilivia.s2.coords.S2Coords.faceSiTitoXYZ
import dilivia.s2.coords.S2Coords.kIJtoPos
import dilivia.s2.coords.S2Coords.kInvertMask
import dilivia.s2.coords.S2Coords.kPosToIJ
import dilivia.s2.coords.S2Coords.kPosToOrientation
import dilivia.s2.coords.S2Coords.kSwapMask
import dilivia.s2.coords.S2Coords.siTiToSt
import dilivia.s2.coords.S2Coords.stToIJ
import dilivia.s2.coords.S2Coords.stToUv
import dilivia.s2.coords.S2Coords.uvToSt
import mu.KotlinLogging
import kotlin.math.max
import kotlin.math.min

// S2PaddedCell represents an S2Cell whose (u,v)-range has been expanded on
// all sides by a given amount of "padding".  Unlike S2Cell, its methods and
// representation are optimized for clipping edges against S2Cell boundaries
// to determine which cells are intersected by a given set of edges.
//
// This class is intended to be copied by value as desired.
class S2PaddedCell {

    val id: S2CellId

    val padding: Double

    // Return the bound for this cell (including padding).
    val bound: MutableR2Rect     // Bound in (u,v)-space

    // The rectangle in (u,v)-space that belongs to all four padded children.
    // It is computed on demand by the middle() accessor method.
    private var middle: R2Rect = R2Rect.empty()

    private val ij_lo: IntArray = IntArray(2)     // Minimum (i,j)-coordinates of this cell, before padding

    val orientation: Int   // Hilbert curve orientation of this cell (see s2coords.h)

    val level: Int;        // Level of this cell (see s2coords.h)

    // Construct an S2PaddedCell for the given cell id and padding.
    constructor(id: S2CellId, padding: Double) {
        this.id = id
        this.padding = padding
        if (id.isFace()) {
            // Fast path for constructing a top-level face (the most common case).
            val limit = 1 + padding
            bound = R2Rect(R1Interval(-limit, limit), R1Interval(-limit, limit)).toMutable()
            middle = R2Rect(R1Interval(-padding, padding), R1Interval(-padding, padding));
            ij_lo[0] = 0
            ij_lo[1] = 0
            orientation = id.face() and 1
            level = 0;
        } else {
            val (_, i, j, orient) = id.toFaceIJOrientation(true)
            level = id.level()
            bound = S2CellId.ijLevelToBoundUV(i, j, level).expanded(padding).toMutable()
            orientation = orient!!
            val ij_size = S2CellId.getSizeIJ(level)
            ij_lo[0] = i and -ij_size
            ij_lo[1] = j and -ij_size
        }
    }

    // Construct the child of "parent" with the given (i,j) index.  The four
    // child cells have indices of (0,0), (0,1), (1,0), (1,1), where the i and j
    // indices correspond to increasing u- and v-values respectively.
    constructor(parent: S2PaddedCell, i: Int, j: Int) {
        logger.trace { "Subdivise cell $parent: i = $i, j = $j " }
        this.padding = parent.padding
        this.bound = parent.bound.toMutable()
        level = parent.level + 1
        // Compute the position and orientation of the child incrementally from the
        // orientation of the parent.
        val pos = kIJtoPos[parent.orientation][2*i+j]
        this.id = parent.id.child(pos)
        val ij_size = S2CellId.getSizeIJ(level)
        ij_lo[0] = parent.ij_lo[0] + i * ij_size
        ij_lo[1] = parent.ij_lo[1] + j * ij_size
        this.orientation = parent.orientation xor kPosToOrientation[pos]
        // For each child, one corner of the bound is taken directly from the parent
        // while the diagonally opposite corner is taken from middle().
        val m = parent.middle()
        logger.trace { "Parent middle: $m" }
        bound[0][1-i] = m[0][1-i]
        bound[1][1-j] = m[1][1-j]
    }

    // Return the "middle" of the padded cell, defined as the rectangle that
    // belongs to all four children.
    //
    // Note that this method is *not* thread-safe, because the return value is
    // computed on demand and cached.  (It is expected that this class will be
    // mainly useful in the context of single-threaded recursive algorithms.)
    fun middle(): R2Rect {
        // We compute this field lazily because it is not needed the majority of the
        // time (i.e., for cells where the recursion terminates).
        if (middle.isEmpty) {
            val ij_size = S2CellId.getSizeIJ(level).toUInt()
            val u = stToUv(siTiToSt((2U * ij_lo[0].toUInt() + ij_size)))
            val v = stToUv(siTiToSt((2U * ij_lo[1].toUInt() + ij_size)))
            middle = R2Rect(R1Interval(u - padding, u + padding), R1Interval(v - padding, v + padding))
        }
        return middle
    }

    // Return the (i,j) coordinates for the child cell at the given traversal
    // position.  The traversal position corresponds to the order in which child
    // cells are visited by the Hilbert curve.
    fun getChildIJ(pos: Int): Pair<Int, Int> {
        val ij = kPosToIJ[orientation][pos]
        val i = ij shr 1
        val j = ij and 1
        return Pair(i, j)
    }

    // Return the smallest cell that contains all descendants of this cell whose
    // bounds intersect "rect".  For algorithms that use recursive subdivision
    // to find the cells that intersect a particular object, this method can be
    // used to skip all the initial subdivision steps where only one child needs
    // to be expanded.
    //
    // Note that this method is not the same as returning the smallest cell that
    // contains the intersection of this cell with "rect".  Because of the
    // padding, even if one child completely contains "rect" it is still
    // possible that a neighboring child also intersects "rect".
    //
    // REQUIRES: bound().Intersects(rect)
    fun shrinkToFit(rect: R2Rect): S2CellId {
        Assertions.assert { bound.intersects(rect) }

        // Quick rejection test: if "rect" contains the center of this cell along
        // either axis, then no further shrinking is possible.
        val ij_size = S2CellId.getSizeIJ(level)
        if (level == 0) {
            // Fast path (most calls to this function start with a face cell).
            if (rect[0].contains(0.0) || rect[1].contains(0.0)) return id
        } else {
            if (rect[0].contains(stToUv(siTiToSt((2 * ij_lo[0] + ij_size).toUInt()))) ||
                    rect[1].contains(stToUv(siTiToSt((2 * ij_lo[1] + ij_size).toUInt())))
            ) {
                return id
            }
        }
        // Otherwise we expand "rect" by the given padding() on all sides and find
        // the range of coordinates that it spans along the i- and j-axes.  We then
        // compute the highest bit position at which the min and max coordinates
        // differ.  This corresponds to the first cell level at which at least two
        // children intersect "rect".

        // Increase the padding to compensate for the error in S2::UVtoST().
        // (The constant below is a provable upper bound on the additional error.)
        val padded = rect.expanded(padding + 1.5 * DBL_EPSILON)
        val ij_min = IntArray(2)  // Min i- or j- coordinate spanned by "padded"
        val ij_xor = IntArray(2)  // XOR of the min and max i- or j-coordinates
        for (d in 0..1) {
            ij_min[d] = max(ij_lo[d], stToIJ(uvToSt(padded[d][0])))
            val ij_max = min(ij_lo[d] + ij_size - 1, stToIJ(uvToSt(padded[d][1])))
            ij_xor[d] = ij_min[d] xor ij_max
        }
        // Compute the highest bit position where the two i- or j-endpoints differ,
        // and then choose the cell level that includes both of these endpoints.  So
        // if both pairs of endpoints are equal we choose kMaxLevel; if they differ
        // only at bit 0, we choose (kMaxLevel - 1), and so on.
        val level_msb = ((ij_xor[0] or ij_xor[1]) shl 1) + 1
        val level = S2CellId.kMaxLevel - Bits.findMSBSetNonZero(level_msb.toUInt())
        if (level <= this.level) return id
        return S2CellId.fromFaceIJ(id.face(), ij_min[0], ij_min[1]).parent(level)
    }

    // Return the center of this cell.
    fun getCenter(): S2Point {
        val ij_size = S2CellId.getSizeIJ(level)
        val si = (2 * ij_lo[0] + ij_size).toUInt()
        val ti = (2 * ij_lo[1] + ij_size).toUInt()
        return faceSiTitoXYZ(id.face(), si, ti).normalize()
    }

    // Return the vertex where the S2 space-filling curve enters this cell.
    fun getEntryVertex(): S2Point {
        // The curve enters at the (0,0) vertex unless the axis directions are
        // reversed, in which case it enters at the (1,1) vertex.
        var i = ij_lo[0].toUInt();
        var j = ij_lo[1].toUInt()
        if ((orientation and kInvertMask) != 0) {
            val ij_size = S2CellId.getSizeIJ(level).toUInt()
            i += ij_size
            j += ij_size
        }
        return faceSiTitoXYZ(id.face(), 2U * i, 2U * j).normalize()
    }

    // Return the vertex where the S2 space-filling curve exits this cell.
    fun getExitVertex(): S2Point {
        // The curve exits at the (1,0) vertex unless the axes are swapped or
        // inverted but not both, in which case it exits at the (0,1) vertex.
        var i = ij_lo[0].toUInt()
        var j = ij_lo[1].toUInt()
        val ij_size = S2CellId.getSizeIJ(level).toUInt()
        if (orientation == 0 || orientation == kSwapMask + kInvertMask) {
            i += ij_size;
        } else {
            j += ij_size;
        }
        return faceSiTitoXYZ(id.face(), 2U * i, 2U * j).normalize()
    }

    override fun toString(): String {
        return "S2PaddedCell(id=$id, padding=$padding)"
    }


    companion object {

        val logger = KotlinLogging.logger {  }

    }

}

