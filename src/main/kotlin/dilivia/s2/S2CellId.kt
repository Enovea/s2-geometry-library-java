/*
 * Copyright 2005 Google Inc.
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

import com.google.common.geometry.MutableInteger
import com.google.common.geometry.S2
import com.google.common.geometry.S2Projections
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.math.R2Point
import java.util.*

/**
 * An S2CellId is a 64-bit unsigned integer that uniquely identifies a
 * cell in the S2 cell decomposition.  It has the following format:
 *
 *   id = [face][face_pos]
 *
 *   face:     a 3-bit number (range 0..5) encoding the cube face.
 *
 *   face_pos: a 61-bit number encoding the position of the center of this
 *             cell along the Hilbert curve over this face (see the Wiki
 *             pages for details).
 *
 * Sequentially increasing cell ids follow a continuous space-filling curve
 * over the entire sphere.  They have the following properties:
 *
 *  - The id of a cell at level k consists of a 3-bit face number followed
 *    by k bit pairs that recursively select one of the four children of
 *    each cell.  The next bit is always 1, and all other bits are 0.
 *    Therefore, the level of a cell is determined by the position of its
 *    lowest-numbered bit that is turned on (for a cell at level k, this
 *    position is 2 * (kMaxLevel - k).)
 *
 *  - The id of a parent cell is at the midpoint of the range of ids spanned
 *    by its children (or by its descendants at any level).
 *
 * Leaf cells are often used to represent points on the unit sphere, and
 * this class provides methods for converting directly between these two
 * representations.  For cells that represent 2D regions rather than
 * discrete point, it is better to use the S2Cell class.
 *
 * @property id The 64-bit unique identifier for this cell.
 * @constructor
 * @param id The id of the cell.
 */
@Strictfp
class S2CellId(val id: Long) : Comparable<S2CellId> {

    // The default constructor returns an invalid cell id.
    constructor() : this(0L)

    // Return the direction vector corresponding to the center of the given
    // cell.  The vector returned by ToPointRaw is not necessarily unit length.
    // This method returns the same result as S2Cell::GetCenter().
    //
    // The maximum directional error in ToPoint() (compared to the exact
    // mathematical result) is 1.5 * DBL_EPSILON radians, and the maximum length
    // error is 2 * DBL_EPSILON (the same as Normalize).
    fun toPoint(): S2Point = normalize(toPointRaw())

    /**
     * Return the direction vector corresponding to the center of the given cell.
     * The vector returned by ToPointRaw is not necessarily unit length.
     */
    fun toPointRaw(): S2Point {
        // First we compute the discrete (i,j) coordinates of a leaf cell contained
        // within the given cell. Given that cells are represented by the Hilbert
        // curve position corresponding at their center, it turns out that the cell
        // returned by ToFaceIJOrientation is always one of two leaf cells closest
        // to the center of the cell (unless the given cell is a leaf cell itself,
        // in which case there is only one possibility).
        //
        // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
        // jmin) be the coordinates of its lower left-hand corner, the leaf cell
        // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
        // (imin + s/2 - 1, jmin + s/2 - 1). We can distinguish these two cases by
        // looking at the low bit of "i" or "j". In the first case the low bit is
        // zero, unless s == 2 (i.e. the level just above leaf cells) in which case
        // the low bit is one.
        //
        // The following calculation converts (i,j) to the (si,ti) coordinates of
        // the cell center. (We need to multiply the coordinates by a factor of 2
        // so that the center of leaf cells can be represented exactly.)
        val (face, i, j, _) = toFaceIJOrientation()
        // System.out.println("i= " + i.intValue() + " j = " + j.intValue());
        val delta = if (isLeaf) 1 else if (i xor (id.toInt() ushr 2) and 1 != 0) 2 else 0
        val si = (i shl 1) + delta - kMaxSize
        val ti = (j shl 1) + delta - kMaxSize
        return faceSiTiToXYZ(face, si, ti)
    }

    /** Return the S2LatLng corresponding to the center of the given cell.  */
    fun toLatLng(): S2LatLng {
        return S2LatLng.fromPoint(toPointRaw())
    }

    // Return the center of the cell in (s,t) coordinates (see s2coords.h).
    fun getCenterST(): R2Point = TODO()

    // Return the edge length of this cell in (s,t)-space.
    fun getSizeST(): Double = getSizeST(level())

    // Return the bounds of this cell in (s,t)-space.
    fun getBoundST(): R2Rect = TODO()

    // Return the center of the cell in (u,v) coordinates (see s2coords.h).
    // Note that the center of the cell is defined as the point at which it is
    // recursively subdivided into four children; in general, it is not at the
    // midpoint of the (u,v) rectangle covered by the cell.
    fun getCenterUV(): R2Point = TODO()

    // Return the bounds of this cell in (u,v)-space.
    fun getBoundUV(): R2Rect = TODO()

    // Return the (face, si, ti) coordinates of the center of the cell.  Note
    // that although (si,ti) coordinates span the range [0,2**31] in general,
    // the cell center coordinates are always in the range [1,2**31-1] and
    // therefore can be represented using a signed 32-bit integer.
    fun getCenterSiTi(): FaceSiTi_TMP {
        // First we compute the discrete (i,j) coordinates of a leaf cell contained
        // within the given cell.  Given that cells are represented by the Hilbert
        // curve position corresponding at their center, it turns out that the cell
        // returned by ToFaceIJOrientation is always one of two leaf cells closest
        // to the center of the cell (unless the given cell is a leaf cell itself,
        // in which case there is only one possibility).
        //
        // Given a cell of size s >= 2 (i.e. not a leaf cell), and letting (imin,
        // jmin) be the coordinates of its lower left-hand corner, the leaf cell
        // returned by ToFaceIJOrientation() is either (imin + s/2, jmin + s/2)
        // (imin + s/2 - 1, jmin + s/2 - 1).  The first case is the one we want.
        // We can distinguish these two cases by looking at the low bit of "i" or
        // "j".  In the second case the low bit is one, unless s == 2 (i.e. the
        // level just above leaf cells) in which case the low bit is zero.
        //
        // In the code below, the expression ((i ^ (int(id_) >> 2)) & 1) is true
        // if we are in the second case described above.
        val (face, i, j, _) = toFaceIJOrientation()
        val delta = if(isLeaf)  1 else if(((i xor (id.toInt() shr 2)) and 1) != 0) 2 else 0;

        // Note that (2 * {i,j} + delta) will never overflow a 32-bit integer.
        return FaceSiTi_TMP(face = face, si = 2 * i + delta, ti = 2 * j + delta)
    }

    // Return true if id() represents a valid cell.
    //
    // All methods require is_valid() to be true unless otherwise specified
    // (although not all methods enforce this).
    val isValid: Boolean
        get() = face() < kNumFaces && lowestOnBit() and 0x1555555555555555L != 0L

    /** Which cube face this cell belongs to, in the range 0..5.  */
    fun face(): Int {
        return (id ushr kPosBits).toInt()
    }

    /**
     * The position of the cell center along the Hilbert curve over this face, in
     * the range 0..(2**kPosBits-1).
     */
    fun pos(): Long {
        return id and (-1L ushr kFaceBits)
    }

    /** Return the subdivision level of the cell (range 0..MAX_LEVEL).  */
    fun level(): Int {
        // Fast path for leaf cells.
        if (isLeaf) {
            return kMaxLevel
        }
        var x = id.toInt()
        var level = -1
        if (x != 0) {
            level += 16
        } else {
            x = (id ushr 32).toInt()
        }
        // We only need to look at even-numbered bits to determine the
        // level of a valid cell id.
        x = x and -x // Get lowest bit.
        if (x and 0x00005555 != 0) {
            level += 8
        }
        if (x and 0x00550055 != 0) {
            level += 4
        }
        if (x and 0x05050505 != 0) {
            level += 2
        }
        if (x and 0x11111111 != 0) {
            level += 1
        }
        // assert (level >= 0 && level <= MAX_LEVEL);
        return level
    }

    // Return the edge length of this cell in (i,j)-space.
    fun getSizeIJ(): Int = getSizeIJ(level())

    /**
     * Return true if this is a leaf cell (more efficient than checking whether
     * level() == MAX_LEVEL).
     */
    val isLeaf: Boolean
        get() = id.toInt() and 1 != 0

    /**
     * Return true if this is a top-level face cell (more efficient than checking
     * whether level() == 0).
     */
    val isFace: Boolean
        get() = id and lowestOnBitForLevel(0) - 1 == 0L

    /**
     * Return the child position (0..3) of this cell's ancestor at the given
     * level, relative to its parent. The argument should be in the range
     * 1..MAX_LEVEL. For example, child_position(1) returns the position of this
     * cell's level-1 ancestor within its top-level face cell.
     */
    fun childPosition(level: Int): Int {
        return (id ushr 2 * (kMaxLevel - level) + 1).toInt() and 3
    }

    // Return the child position (0..3) of this cell within its parent.
    // REQUIRES: level() >= 1.
    fun childPosition(): Int = TODO()

    // These methods return the range of cell ids that are contained within this
    // cell (including itself).  The range is *inclusive* (i.e. test using >=
    // and <=) and the return values of both methods are valid leaf cell ids.
    // In other words, a.contains(b) if and only if
    //
    //     (b >= a.range_min() && b <= a.range_max())
    //
    // If you want to iterate through all the descendants of this cell at a
    // particular level, use child_begin(level) and child_end(level) instead.
    // Also see maximum_tile(), which can be used to iterate through a range of
    // cells using S2CellIds at different levels that are as large as possible.
    //
    // If you need to convert the range to a semi-open interval [min, limit)
    // (e.g., in order to use a key-value store that only supports semi-open
    // range queries), do not attempt to define "limit" as range_max.next().
    // The problem is that leaf S2CellIds are 2 units apart, so the semi-open
    // interval [min, limit) includes an additional value (range_max.id() + 1)
    // which is happens to be a valid S2CellId about one-third of the time and
    // is *never* contained by this cell.  (It always correpsonds to a cell that
    // is larger than this one.)  You can define "limit" as (range_max.id() + 1)
    // if necessary (which is not always a valid S2CellId but can still be used
    // with FromToken/ToToken), or you can convert range_max() to the key space
    // of your key-value store and define "limit" as Successor(key).
    //
    // Note that Sentinel().range_min() == Sentinel.range_max() == Sentinel().

    fun rangeMin(): S2CellId {
        return S2CellId(id - (lowestOnBit() - 1))
    }

    fun rangeMax(): S2CellId {
        return S2CellId(id + (lowestOnBit() - 1))
    }

    /** Return true if the given cell is contained within this one.  */
    operator fun contains(other: S2CellId): Boolean {
        // assert (isValid() && other.isValid());
        return other.greaterOrEquals(rangeMin()) && other.lessOrEquals(rangeMax())
    }

    /** Return true if the given cell intersects this one.  */
    fun intersects(other: S2CellId): Boolean {
        // assert (isValid() && other.isValid());
        return (other.rangeMin().lessOrEquals(rangeMax())
                && other.rangeMax().greaterOrEquals(rangeMin()))
    }

    // Return the cell at the previous level or at the given level (which must
    // be less than or equal to the current level).
    fun parent(): S2CellId {
        // assert (isValid() && level() > 0);
        val newLsb = lowestOnBit() shl 2
        return S2CellId(id and -newLsb or newLsb)
    }

    fun parent(level: Int): S2CellId {
        // assert (isValid() && level >= 0 && level <= this.level());
        val newLsb = lowestOnBitForLevel(level)
        return S2CellId(id and -newLsb or newLsb)
    }

    // Return the immediate child of this cell at the given traversal order
    // position (in the range 0 to 3).  This cell must not be a leaf cell.
    fun child(position: Int): S2CellId = TODO()

    // Iterator-style methods for traversing the immediate children of a cell or
    // all of the children at a given level (greater than or equal to the current
    // level).  Note that the end value is exclusive, just like standard STL
    // iterators, and may not even be a valid cell id.  You should iterate using
    // code like this:
    //
    //   for(S2CellId c = id.child_begin(); c != id.child_end(); c = c.next())
    //     ...
    //
    // The convention for advancing the iterator is "c = c.next()" rather
    // than "++c" to avoid possible confusion with incrementing the
    // underlying 64-bit cell id.
    fun childBegin(): S2CellId {
        // assert (isValid() && level() < MAX_LEVEL);
        val oldLsb = lowestOnBit()
        return S2CellId(id - oldLsb + (oldLsb ushr 2))
    }

    fun childBegin(level: Int): S2CellId {
        // assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
        return S2CellId(id - lowestOnBit() + lowestOnBitForLevel(level))
    }

    fun childEnd(): S2CellId {
        // assert (isValid() && level() < MAX_LEVEL);
        val oldLsb = lowestOnBit()
        return S2CellId(id + oldLsb + (oldLsb ushr 2))
    }

    fun childEnd(level: Int): S2CellId {
        // assert (isValid() && level >= this.level() && level <= MAX_LEVEL);
        return S2CellId(id + lowestOnBit() + lowestOnBitForLevel(level))
    }

    // Return the next/previous cell at the same level along the Hilbert curve.
    // Works correctly when advancing from one face to the next, but
    // does *not* wrap around from the last face to the first or vice versa.
    /**
     * Return the next cell at the same level along the Hilbert curve. Works
     * correctly when advancing from one face to the next, but does *not* wrap
     * around from the last face to the first or vice versa.
     */
    operator fun next(): S2CellId {
        return S2CellId(id + (lowestOnBit() shl 1))
    }

    /**
     * Return the previous cell at the same level along the Hilbert curve. Works
     * correctly when advancing from one face to the next, but does *not* wrap
     * around from the last face to the first or vice versa.
     */
    fun prev(): S2CellId {
        return S2CellId(id - (lowestOnBit() shl 1))
    }

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position is never advanced past End() or before Begin().
    fun advance(steps: Long): S2CellId {
        if (steps == 0L) return this

        // We clamp the number of steps if necessary to ensure that we do not
        // advance past the End() or before the Begin() of this level.  Note that
        // min_steps and max_steps always fit in a signed 64-bit integer.
        var s = steps
        val step_shift = 2 * (kMaxLevel - level()) + 1
        if (s < 0) {
            val min_steps = -(id shr step_shift)
            if (s < min_steps) s = min_steps
        } else {
            val max_steps = (kWrapOffset + lowestOnBit() - id) shr step_shift
            if (s > max_steps) s = max_steps
        }
        // If steps is negative, then shifting it left has undefined behavior.
        // Cast to uint64 for a 2's complement answer.
        return S2CellId((id.toULong() + (s.toULong() shl step_shift)).toLong())
    }

    // Returns the number of steps that this cell is from Begin(level()). The
    // return value is always non-negative.
    fun distanceFromBegin(): Long = TODO()

    /**
     * Like next(), but wraps around from the last face to the first and vice
     * versa. Should *not* be used for iteration in conjunction with
     * child_begin(), child_end(), Begin(), or End().
     */
    fun nextWrap(): S2CellId {
        val n = next()
        return if (unsignedLongLessThan(n.id, kWrapOffset)) {
            n
        } else S2CellId(n.id - kWrapOffset)
    }

    /**
     * Like prev(), but wraps around from the last face to the first and vice
     * versa. Should *not* be used for iteration in conjunction with
     * child_begin(), child_end(), Begin(), or End().
     */
    fun prevWrap(): S2CellId {
        val p = prev()
        return if (p.id < kWrapOffset) {
            p
        } else S2CellId(p.id + kWrapOffset)
    }

    // This method advances or retreats the indicated number of steps along the
    // Hilbert curve at the current level, and returns the new position.  The
    // position wraps between the first and last faces as necessary.  The input
    // must be a valid cell id.
    fun advanceWrap(steps: Long): S2CellId = TODO()

    // Return the largest cell with the same range_min() and such that
    // range_max() < limit.range_min().  Returns "limit" if no such cell exists.
    // This method can be used to generate a small set of S2CellIds that covers
    // a given range (a "tiling").  This example shows how to generate a tiling
    // for a semi-open range of leaf cells [start, limit):
    //
    //   for (S2CellId id = start.maximum_tile(limit);
    //        id != limit; id = id.next().maximum_tile(limit)) { ... }
    //
    // Note that in general the cells in the tiling will be of different sizes;
    // they gradually get larger (near the middle of the range) and then
    // gradually get smaller (as "limit" is approached).
    fun maximumTile(limit: S2CellId): S2CellId = TODO()

    // Methods to encode and decode cell ids to compact text strings suitable
    // for display or indexing.  Cells at lower levels (i.e. larger cells) are
    // encoded into fewer characters.  The maximum token length is 16.
    //
    // Tokens preserve ordering, i.e. ToToken(x) < ToToken(y) iff x < y.
    //
    // ToToken() returns a string by value for convenience; the compiler
    // does this without intermediate copying in most cases.
    //
    // These methods guarantee that FromToken(ToToken(x)) == x even when
    // "x" is an invalid cell id.  All tokens are alphanumeric strings.
    // FromToken() returns S2CellId::None() for malformed inputs.
    fun toToken(): String {
        if (id == 0L) {
            return "X"
        }
        val hex = java.lang.Long.toHexString(id).toLowerCase(Locale.ENGLISH)
        val sb = StringBuilder(16)
        for (i in hex.length..15) {
            sb.append('0')
        }
        sb.append(hex)
        for (len in 16 downTo 1) {
            if (sb[len - 1] != '0') {
                return sb.substring(0, len)
            }
        }
        throw RuntimeException("Shouldn't make it here")
    }

    // Creates a human readable debug string.  Used for << and available for
    // direct usage as well.  The format is "f/dd..d" where "f" is a digit in
    // the range [0-5] representing the S2CellId face, and "dd..d" is a string
    // of digits in the range [0-3] representing each child's position with
    // respect to its parent.  (Note that the latter string may be empty.)
    //
    // For example "4/" represents S2CellId::FromFace(4), and "3/02" represents
    // S2CellId::FromFace(3).child(0).child(2).
    override fun toString(): String {
        return ("(face=" + face() + ", pos=" + java.lang.Long.toHexString(pos()) + ", level=" + level() + ")")
    }

    /**
     * Return the four cells that are adjacent across the cell's four edges.
     * Neighbors are returned in the order defined by S2Cell::GetEdge. All
     * neighbors are guaranteed to be distinct.
     */
    fun getEdgeNeighbors(neighbors: Array<S2CellId>) {
        val level = level()
        val size = 1 shl kMaxLevel - level
        val (face, i, j, _) = toFaceIJOrientation()

        // Edges 0, 1, 2, 3 are in the S, E, N, W directions.
        neighbors[0] = fromFaceIJSame(face, i, j - size, j - size >= 0).parent(level)
        neighbors[1] = fromFaceIJSame(face, i + size, j, i + size < kMaxSize).parent(level)
        neighbors[2] = fromFaceIJSame(face, i, j + size, j + size < kMaxSize).parent(level)
        neighbors[3] = fromFaceIJSame(face, i - size, j, i - size >= 0).parent(level)
    }

    /**
     * Return the neighbors of closest vertex to this cell at the given level, by
     * appending them to "output". Normally there are four neighbors, but the
     * closest vertex may only have three neighbors if it is one of the 8 cube
     * vertices.
     *
     * Requires: level < this.evel(), so that we can determine which vertex is
     * closest (in particular, level == MAX_LEVEL is not allowed).
     */
    fun appendVertexNeighbors(level: Int, output: MutableList<S2CellId>) {
        // "level" must be strictly less than this cell's level so that we can
        // determine which vertex this cell is closest to.
        // assert (level < this.level());
        val (face, i, j, _) = toFaceIJOrientation()

        // Determine the i- and j-offsets to the closest neighboring cell in each
        // direction. This involves looking at the next bit of "i" and "j" to
        // determine which quadrant of this->parent(level) this cell lies in.
        val halfsize = 1 shl kMaxLevel - (level + 1)
        val size = halfsize shl 1
        val isame: Boolean
        val jsame: Boolean
        val ioffset: Int
        val joffset: Int
        if (i and halfsize != 0) {
            ioffset = size
            isame = i + size < kMaxSize
        } else {
            ioffset = -size
            isame = i - size >= 0
        }
        if (j and halfsize != 0) {
            joffset = size
            jsame = j + size < kMaxSize
        } else {
            joffset = -size
            jsame = j - size >= 0
        }
        output.add(parent(level))
        output.add(fromFaceIJSame(face, i + ioffset, j, isame).parent(level))
        output.add(fromFaceIJSame(face, i, j + joffset, jsame).parent(level))
        // If i- and j- edge neighbors are *both* on a different face, then this
        // vertex only has three neighbors (it is one of the 8 cube vertices).
        if (isame || jsame) {
            output.add(fromFaceIJSame(face, i + ioffset, j + joffset, isame && jsame).parent(level))
        }
    }

    /**
     * Append all neighbors of this cell at the given level to "output". Two cells
     * X and Y are neighbors if their boundaries intersect but their interiors do
     * not. In particular, two cells that intersect at a single point are
     * neighbors.
     *
     * Requires: nbr_level >= this->level(). Note that for cells adjacent to a
     * face vertex, the same neighbor may be appended more than once.
     */
    fun appendAllNeighbors(nbrLevel: Int, output: MutableList<S2CellId>) {
        var (face, i, j, _) = toFaceIJOrientation(false)

        // Find the coordinates of the lower left-hand leaf cell. We need to
        // normalize (i,j) to a known position within the cell because nbr_level
        // may be larger than this cell's level.
        val size = 1 shl kMaxLevel - level()
        i = (i and -size)
        j = (j and -size)
        val nbrSize = 1 shl kMaxLevel - nbrLevel
        // assert (nbrSize <= size);

        // We compute the N-S, E-W, and diagonal neighbors in one pass.
        // The loop test is at the end of the loop to avoid 32-bit overflow.
        var k = -nbrSize
        while (true) {
            var sameFace: Boolean
            if (k < 0) {
                sameFace = j + k >= 0
            } else if (k >= size) {
                sameFace = j + k < kMaxSize
            } else {
                sameFace = true
                // North and South neighbors.
                output.add(fromFaceIJSame(face, i + k, j - nbrSize, j - size >= 0).parent(nbrLevel))
                output.add(fromFaceIJSame(face, i + k, j + size, j + size < kMaxSize).parent(nbrLevel))
            }
            // East, West, and Diagonal neighbors.
            output.add(fromFaceIJSame(face, i - nbrSize, j + k, sameFace && i - size >= 0).parent(nbrLevel))
            output.add(fromFaceIJSame(face, i + size, j + k, sameFace && i + size < kMaxSize).parent(nbrLevel))
            if (k >= size) {
                break
            }
            k += nbrSize
        }
    }

    /**
     * Return the (face, i, j) coordinates for the leaf cell corresponding to this
     * cell id. Since cells are represented by the Hilbert curve position at the
     * center of the cell, the returned (i,j) for non-leaf cells will be a leaf
     * cell adjacent to the cell center. If "orientation" is non-NULL, also return
     * the Hilbert curve orientation for the current cell.
     */
    @JvmOverloads
    fun toFaceIJOrientation(computeOrientation: Boolean = false): FaceIJ {
        // System.out.println("Entering toFaceIjorientation");
        val face = face()
        var bits = face and LookupCellTables.kSwapMask
        val pi = MutableInteger(0)
        val pj = MutableInteger(0)

        // System.out.println("face = " + face + " bits = " + bits);

        // Each iteration maps 8 bits of the Hilbert curve position into
        // 4 bits of "i" and "j". The lookup table transforms a key of the
        // form "ppppppppoo" to a value of the form "iiiijjjjoo", where the
        // letters [ijpo] represents bits of "i", "j", the Hilbert curve
        // position, and the Hilbert curve orientation respectively.
        //
        // On the first iteration we need to be careful to clear out the bits
        // representing the cube face.
        for (k in 7 downTo 0) {
            bits = getBits1(pi, pj, k, bits)
            // System.out.println("pi = " + pi + " pj= " + pj + " bits = " + bits);
        }
        var orientation: Int? = null
        if (computeOrientation) {
            // The position of a non-leaf cell at level "n" consists of a prefix of
            // 2*n bits that identifies the cell, followed by a suffix of
            // 2*(MAX_LEVEL-n)+1 bits of the form 10*. If n==MAX_LEVEL, the suffix is
            // just "1" and has no effect. Otherwise, it consists of "10", followed
            // by (MAX_LEVEL-n-1) repetitions of "00", followed by "0". The "10" has
            // no effect, while each occurrence of "00" has the effect of reversing
            // the kSwapMask bit.
            // assert (S2.POS_TO_ORIENTATION[2] == 0);
            // assert (S2.POS_TO_ORIENTATION[0] == S2.SWAP_MASK);
            if (lowestOnBit() and 0x1111111111111110L != 0L) {
                bits = bits xor S2.SWAP_MASK
            }
            orientation = bits
        }
        return FaceIJ(face, pi.intValue(), pj.intValue(), orientation)
    }

    // Return the lowest-numbered bit that is on for this cell id, which is
    // equal to (uint64{1} << (2 * (kMaxLevel - level))).  So for example,
    // a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the
    // first test is more efficient.
    /** Return the lowest-numbered bit that is on for cells at the given level.  */
    fun lowestOnBit(): Long {
        return id and -id
    }

    private fun getBits1(i: MutableInteger, j: MutableInteger, k: Int, bits: Int): Int {
        var bits = bits
        val nbits = if (k == 7) kMaxLevel - 7 * LookupCellTables.kLookupBits else LookupCellTables.kLookupBits
        bits += (id ushr k * 2 * LookupCellTables.kLookupBits + 1).toInt() and
                (1 shl 2 * nbits) - 1 shl 2
        /*
     * System.out.println("id is: " + id_); System.out.println("bits is " +
     * bits); System.out.println("lookup_ij[bits] is " + lookup_ij[bits]);
     */bits = LookupCellTables.lookupIj[bits]
        i.setValue(i.intValue()
                + (bits shr LookupCellTables.kLookupBits + 2 shl k * LookupCellTables.kLookupBits))
        /*
     * System.out.println("left is " + ((bits >> 2) & ((1 << kLookupBits) -
     * 1))); System.out.println("right is " + (k * kLookupBits));
     * System.out.println("j is: " + j.intValue()); System.out.println("addition
     * is: " + ((((bits >> 2) & ((1 << kLookupBits) - 1))) << (k *
     * kLookupBits)));
     */j.setValue(j.intValue()
                + (bits shr 2 and (1 shl LookupCellTables.kLookupBits) - 1 shl k * LookupCellTables.kLookupBits))
        bits = bits and (LookupCellTables.kSwapMask or LookupCellTables.kInvertMask)
        return bits
    }


    override fun equals(that: Any?): Boolean {
        if (that !is S2CellId) {
            return false
        }
        return id == that.id
    }

    override fun hashCode(): Int {
        return ((id ushr 32) + id).toInt()
    }

    fun lessThan(x: S2CellId): Boolean {
        return unsignedLongLessThan(id, x.id)
    }

    fun greaterThan(x: S2CellId): Boolean {
        return unsignedLongGreaterThan(id, x.id)
    }

    fun lessOrEquals(x: S2CellId): Boolean {
        return unsignedLongLessThan(id, x.id) || id == x.id
    }

    fun greaterOrEquals(x: S2CellId): Boolean {
        return unsignedLongGreaterThan(id, x.id) || id == x.id
    }



    override fun compareTo(that: S2CellId): Int {
        return if (unsignedLongLessThan(id, that.id)) -1 else if (unsignedLongGreaterThan(id, that.id)) 1 else 0
    }

    companion object {

        // The extra position bit (61 rather than 60) let us encode each cell as its
        // Hilbert curve position at the cell center (which is halfway along the
        // portion of the Hilbert curve that fills that cell).
        const val kFaceBits = 3
        const val kNumFaces = 6
        const val kMaxLevel = 30 // Valid levels: 0..kMaxLevel
        const val kPosBits = 2 * kMaxLevel + 1
        const val kMaxSize = 1 shl kMaxLevel

        // Constant related to unsigned long's
        const val kMaxUnsigned = -1L // Equivalent to 0xffffffffffffffffL

        /** The default constructor returns an invalid cell id.  */
        @JvmStatic
        fun none(): S2CellId = S2CellId()

        /**
         * Returns an invalid cell id guaranteed to be larger than any valid cell id.
         * Useful for creating indexes.
         */
        @JvmStatic
        fun sentinel(): S2CellId = S2CellId(kMaxUnsigned) // -1

        // Return the cell corresponding to a given S2 cube face.
        @JvmStatic
        fun fromFace(face: Int): S2CellId = S2CellId((face.toLong() shl kPosBits) + lowestOnBitForLevel(0))

        /**
         * Return a cell given its face (range 0..5), 61-bit Hilbert curve position
         * within that face, and level (range 0..MAX_LEVEL). The given position will
         * be modified to correspond to the Hilbert curve position at the center of
         * the returned cell. This is a static function rather than a constructor in
         * order to give names to the arguments.
         */
        @JvmStatic
        fun fromFacePosLevel(face: Int, pos: Long, level: Int): S2CellId = S2CellId((face.toLong() shl kPosBits) + (pos or 1)).parent(level)

        // Return the edge length in (s,t)-space of cells at the given level.
        fun getSizeST(level: Int): Double = TODO()

        // Expand a rectangle in (u,v)-space so that it contains all points within
        // the given distance of the boundary, and return the smallest such
        // rectangle.  If the distance is negative, then instead shrink this
        // rectangle so that it excludes all points within the given absolute
        // distance of the boundary.
        //
        // Distances are measured *on the sphere*, not in (u,v)-space.  For example,
        // you can use this method to expand the (u,v)-bound of an S2CellId so that
        // it contains all points within 5km of the original cell.  You can then
        // test whether a point lies within the expanded bounds like this:
        //
        //   R2Point uv;
        //   if (S2::FaceXYZtoUV(face, point, &uv) && bound.Contains(uv)) { ... }
        //
        // Limitations:
        //
        //  - Because the rectangle is drawn on one of the six cube-face planes
        //    (i.e., {x,y,z} = +/-1), it can cover at most one hemisphere.  This
        //    limits the maximum amount that a rectangle can be expanded.  For
        //    example, S2CellId bounds can be expanded safely by at most 45 degrees
        //    (about 5000 km on the Earth's surface).
        //
        //  - The implementation is not exact for negative distances.  The resulting
        //    rectangle will exclude all points within the given distance of the
        //    boundary but may be slightly smaller than necessary.
        fun expandedByDistanceUV(uv: R2Rect, distance: S1Angle): R2Rect = TODO()

        // Like the above, but return the size of cells at the given level.
        fun getSizeIJ(level: Int): Int = 1 shl (kMaxLevel - level)

        // Iterator-style methods for traversing all the cells along the Hilbert
        // curve at a given level (across all 6 faces of the cube).  Note that the
        // end value is exclusive (just like standard STL iterators), and is not a
        // valid cell id.
        @JvmStatic
        fun begin(level: Int): S2CellId {
            return fromFacePosLevel(0, 0, 0).childBegin(level)
        }

        @JvmStatic
        fun end(level: Int): S2CellId {
            return fromFacePosLevel(5, 0, 0).childEnd(level)
        }

        /**
         * Decodes the cell id from a compact text string suitable for display or
         * indexing. Cells at lower levels (i.e. larger cells) are encoded into
         * fewer characters. The maximum token length is 16.
         *
         * @param token the token to decode
         * @return the S2CellId for that token
         * @throws NumberFormatException if the token is not formatted correctly
         */
        @JvmStatic
        fun fromToken(token: String): S2CellId {
            if (token.length == 0) {
                throw NumberFormatException("Empty string in S2CellId.fromToken")
            }
            if (token.length > 16 || "X" == token) {
                return none()
            }
            var value: Long = 0
            for (pos in 0..15) {
                var digit = 0
                if (pos < token.length) {
                    digit = Character.digit(token[pos], 16)
                    if (digit == -1) {
                        throw NumberFormatException(token)
                    }
                    if (overflowInParse(value, digit)) {
                        throw NumberFormatException("Too large for unsigned long: $token")
                    }
                }
                value = value * 16 + digit
            }
            return S2CellId(value)
        }

        // Converts a string in the format returned by ToString() to an S2CellId.
        // Returns S2CellId::None() if the string could not be parsed.
        //
        // The method name includes "Debug" in order to avoid possible confusion
        // with FromToken() above.
        fun fromDebugString(str: String): S2CellId = TODO()

        // ///////////////////////////////////////////////////////////////////
        // Low-level methods.
        /**
         * Return a leaf cell given its cube face (range 0..5) and i- and
         * j-coordinates (see s2.h).
         */
        @JvmStatic
        fun fromFaceIJ(faceIJ: FaceIJ): S2CellId = fromFaceIJ(faceIJ.face, faceIJ.i, faceIJ.j)

        @JvmStatic
        fun fromFaceIJ(face: Int, i: Int, j: Int): S2CellId {
            // Optimization notes:
            // - Non-overlapping bit fields can be combined with either "+" or "|".
            // Generally "+" seems to produce better code, but not always.

            // gcc doesn't have very good code generation for 64-bit operations.
            // We optimize this by computing the result as two 32-bit integers
            // and combining them at the end. Declaring the result as an array
            // rather than local variables helps the compiler to do a better job
            // of register allocation as well. Note that the two 32-bits halves
            // get shifted one bit to the left when they are combined.
            val n = longArrayOf(0, (face shl (kPosBits - 33).toLong().toInt()).toLong())

            // Alternating faces have opposite Hilbert curve orientations; this
            // is necessary in order for all faces to have a right-handed
            // coordinate system.
            var bits = face and LookupCellTables.kSwapMask

            // Each iteration maps 4 bits of "i" and "j" into 8 bits of the Hilbert
            // curve position. The lookup table transforms a 10-bit key of the form
            // "iiiijjjjoo" to a 10-bit value of the form "ppppppppoo", where the
            // letters [ijpo] denote bits of "i", "j", Hilbert curve position, and
            // Hilbert curve orientation respectively.
            for (k in 7 downTo 0) {
                bits = LookupCellTables.getBits(n, i, j, k, bits)
            }
            return S2CellId(((n[1] shl 32) + n[0] shl 1) + 1)
        }

        /**
         * Return the lowest-numbered bit that is on for this cell id, which is equal
         * to (uint64(1) << (2 * (MAX_LEVEL - level))). So for example, a.lsb() <=
         * b.lsb() if and only if a.level() >= b.level(), but the first test is more
         * efficient.
         */
        @JvmStatic
        fun lowestOnBitForLevel(level: Int): Long {
            return 1L shl (2 * (kMaxLevel - level))
        }

        // Return the bound in (u,v)-space for the cell at the given level containing
        // the leaf cell with the given (i,j)-coordinates.
        fun ijLevelToBoundUV(ij: IntArray, level: Int): R2Rect = TODO()

        /**
         * This is the offset required to wrap around from the beginning of the
         * Hilbert curve to the end or vice versa; see next_wrap() and prev_wrap().
         */
        private const val kWrapOffset = kNumFaces.toLong() shl kPosBits


        // Given a face and a point (i,j) where either i or j is outside the valid
        // range [0..kMaxSize-1], this function first determines which neighboring
        // face "contains" (i,j), and then returns the leaf cell on that face which
        // is adjacent to the given face and whose distance from (i,j) is minimal.
        private fun fromFaceIJWrap(face: Int, i: Int, j: Int): S2CellId {
            // Convert i and j to the coordinates of a leaf cell just beyond the
            // boundary of this face. This prevents 32-bit overflow in the case
            // of finding the neighbors of a face cell, and also means that we
            // don't need to worry about the distinction between (s,t) and (u,v).
            var face = face
            var i = i
            var j = j
            i = Math.max(-1, Math.min(kMaxSize, i))
            j = Math.max(-1, Math.min(kMaxSize, j))

            // Find the (s,t) coordinates corresponding to (i,j). At least one
            // of these coordinates will be just outside the range [0, 1].
            val kScale = 1.0 / kMaxSize
            val s = kScale * ((i shl 1) + 1 - kMaxSize)
            val t = kScale * ((j shl 1) + 1 - kMaxSize)

            // Find the leaf cell coordinates on the adjacent face, and convert
            // them to a cell id at the appropriate level.
            val p = S2Projections.faceUvToXyz(face, s, t)
            face = S2Projections.xyzToFace(p)
            val st = S2Projections.validFaceXyzToUv(face, p)
            return fromFaceIJ(face, stToIJ(st.x()), stToIJ(st.y()))
        }

        /**
         * Public helper function that calls FromFaceIJ if sameFace is true, or
         * FromFaceIJWrap if sameFace is false.
         */
        fun fromFaceIJSame(face: Int, i: Int, j: Int, sameFace: Boolean): S2CellId {
            return if (sameFace) {
                fromFaceIJ(face, i, j)
            } else {
                fromFaceIJWrap(face, i, j)
            }
        }



        /**
         * Construct a leaf cell containing the given point "p".  Usually there is
         * is exactly one such cell, but for points along the edge of a cell, any
         * adjacent cell may be (deterministically) chosen.  This is because
         * S2CellIds are considered to be closed sets.  The returned cell will
         * always contain the given point, i.e.
         *
         *   S2Cell(S2CellId(p)).Contains(p)
         *
         * is always true.  The point "p" does not need to be normalized.
         *
         * If instead you want every point to be contained by exactly one S2Cell,
         * you will need to convert the S2CellIds to S2Loops (which implement point
         * containment this way).
         */
        @JvmStatic
        fun fromPoint(p: S2Point): S2CellId {
            val face = S2Projections.xyzToFace(p)
            val uv = S2Projections.validFaceXyzToUv(face, p)
            val i = stToIJ(S2Projections.uvToST(uv.x()))
            val j = stToIJ(S2Projections.uvToST(uv.y()))
            return fromFaceIJ(face, i, j)
        }

        /** Return the leaf cell containing the given S2LatLng.  */
        @JvmStatic
        fun fromLatLng(ll: S2LatLng): S2CellId {
            return fromPoint(ll.toPoint())
        }


        /**
         * Returns true if (current * radix) + digit is a number too large to be
         * represented by an unsigned long.  This is useful for detecting overflow
         * while parsing a string representation of a number.
         * Does not verify whether supplied radix is valid, passing an invalid radix
         * will give undefined results or an ArrayIndexOutOfBoundsException.
         */
        /**
         * Returns true if (current * 10) + digit is a number too large to be
         * represented by an unsigned long.  This is useful for detecting overflow
         * while parsing a string representation of a number.
         */
        private fun overflowInParse(current: Long, digit: Int, radix: Int = 10): Boolean {
            if (current >= 0) {
                if (current < maxValueDivs[radix]) {
                    return false
                }
                return if (current > maxValueDivs[radix]) {
                    true
                } else digit > maxValueMods[radix]
                // current == maxValueDivs[radix]
            }

            // current < 0: high bit is set
            return true
        }

        // calculated as 0xffffffffffffffff / radix
        private val maxValueDivs = longArrayOf(0, 0,  // 0 and 1 are invalid
                9223372036854775807L, 6148914691236517205L, 4611686018427387903L,  // 2-4
                3689348814741910323L, 3074457345618258602L, 2635249153387078802L,  // 5-7
                2305843009213693951L, 2049638230412172401L, 1844674407370955161L,  // 8-10
                1676976733973595601L, 1537228672809129301L, 1418980313362273201L,  // 11-13
                1317624576693539401L, 1229782938247303441L, 1152921504606846975L,  // 14-16
                1085102592571150095L, 1024819115206086200L, 970881267037344821L,  // 17-19
                922337203685477580L, 878416384462359600L, 838488366986797800L,  // 20-22
                802032351030850070L, 768614336404564650L, 737869762948382064L,  // 23-25
                709490156681136600L, 683212743470724133L, 658812288346769700L,  // 26-28
                636094623231363848L, 614891469123651720L, 595056260442243600L,  // 29-31
                576460752303423487L, 558992244657865200L, 542551296285575047L,  // 32-34
                527049830677415760L, 512409557603043100L) // 35-36

        // calculated as 0xffffffffffffffff % radix
        private val maxValueMods = intArrayOf(0, 0,  // 0 and 1 are invalid
                1, 0, 3, 0, 3, 1, 7, 6, 5, 4, 3, 2, 1, 0, 15, 0, 15, 16, 15, 15,  // 2-21
                15, 5, 15, 15, 15, 24, 15, 23, 15, 15, 31, 15, 17, 15, 15) // 22-36


        /**
         * Return the i- or j-index of the leaf cell containing the given s- or
         * t-value.
         */
        private fun stToIJ(s: Double): Int {
            // Converting from floating-point to integers via static_cast is very slow
            // on Intel processors because it requires changing the rounding mode.
            // Rounding to the nearest integer using FastIntRound() is much faster.
            val m = kMaxSize / 2 // scaling multiplier
            return Math
                    .max(0, Math.min(2 * m - 1.toLong(), Math.round(m * s + (m - 0.5)))).toInt()
        }

        /**
         * Convert (face, si, ti) coordinates (see s2.h) to a direction vector (not
         * necessarily unit length).
         */
        private fun faceSiTiToXYZ(face: Int, si: Int, ti: Int): S2Point {
            val kScale = 1.0 / kMaxSize
            val u = S2Projections.stToUV(kScale * si)
            val v = S2Projections.stToUV(kScale * ti)
            return S2Projections.faceUvToXyz(face, u, v)
        }

        /**
         * Returns true if x1 < x2, when both values are treated as unsigned.
         */
        fun unsignedLongLessThan(x1: Long, x2: Long): Boolean {
            return x1 + Long.MIN_VALUE < x2 + Long.MIN_VALUE
        }

        /**
         * Returns true if x1 > x2, when both values are treated as unsigned.
         */
        fun unsignedLongGreaterThan(x1: Long, x2: Long): Boolean {
            return x1 + Long.MIN_VALUE > x2 + Long.MIN_VALUE
        }


    }

}


data class FaceSiTi_TMP(val face: Int, val si: Int, val ti: Int)

data class FaceIJ(val face: Int, val i: Int, val j: Int, val orientation: Int?)

// The following lookup tables are used to convert efficiently between an
// (i,j) cell index and the corresponding position along the Hilbert curve.
// "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
// orientation of the current cell into 8 bits representing the order in which
// that subcell is visited by the Hilbert curve, plus 2 bits indicating the
// new orientation of the Hilbert curve within that subcell. (Cell
// orientations are represented as combination of kSwapMask and kInvertMask.)
//
// "lookup_ij" is an inverted table used for mapping in the opposite
// direction.
//
// We also experimented with looking up 16 bits at a time (14 bits of position
// plus 2 of orientation) but found that smaller lookup tables gave better
// performance. (2KB fits easily in the primary cache.)
// Values for these constants are *declared* in the *.h file. Even though
// the declaration specifies a value for the constant, that declaration
// is not a *definition* of storage for the value. Because the values are
// supplied in the declaration, we don't need the values here. Failing to
// define storage causes link errors for any code that tries to take the
// address of one of these values.
object LookupCellTables {

    const val kLookupBits = 4
    const val kSwapMask = 0x01
    const val kInvertMask = 0x02
    val lookupPos = IntArray(1 shl 2 * kLookupBits + 2)
    val lookupIj = IntArray(1 shl 2 * kLookupBits + 2)

    init {
        initLookupCell(0, 0, 0, 0, 0, 0)
        initLookupCell(0, 0, 0, kSwapMask, 0, kSwapMask)
        initLookupCell(0, 0, 0, kInvertMask, 0, kInvertMask)
        initLookupCell(0, 0, 0, kSwapMask or kInvertMask, 0, kSwapMask or kInvertMask)
    }

    private fun initLookupCell(level: Int, i: Int, j: Int, origOrientation: Int, pos: Int, orientation: Int) {
        var currentLevel = level
        var currentI = i
        var currentJ = j
        var currentPos = pos
        if (currentLevel == kLookupBits) {
            val ij = (currentI shl kLookupBits) + currentJ
            lookupPos[(ij shl 2) + origOrientation] = (currentPos shl 2) + orientation
            lookupIj[(currentPos shl 2) + origOrientation] = (ij shl 2) + orientation
        } else {
            currentLevel++
            currentI = currentI shl 1
            currentJ = currentJ shl 1
            currentPos = currentPos shl 2
            // Initialize each sub-cell recursively.
            for (subPos in 0..3) {
                val ij = S2.posToIJ(orientation, subPos)
                val orientationMask = S2.posToOrientation(subPos)
                initLookupCell(currentLevel, currentI + (ij ushr 1), currentJ + (ij and 1), origOrientation,
                        currentPos + subPos, orientation xor orientationMask)
            }
        }
    }

    fun getBits(n: LongArray, i: Int, j: Int, k: Int, bits: Int): Int {
        var bits = bits
        val mask = (1 shl kLookupBits) - 1
        bits += i shr k * kLookupBits and mask shl kLookupBits + 2
        bits += j shr k * kLookupBits and mask shl 2
        bits = lookupPos[bits]
        n[k shr 2] = n[k shr 2] or (bits.toLong() shr 2 shl (k and 3) * 2 * kLookupBits)
        bits = bits and (kSwapMask or kInvertMask)
        return bits
    }

}