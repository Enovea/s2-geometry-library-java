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

/**
 * The following lookup tables are used to convert efficiently between an
 * (i,j) cell index and the corresponding position along the Hilbert curve.
 * "lookup_pos" maps 4 bits of "i", 4 bits of "j", and 2 bits representing the
 * orientation of the current cell into 8 bits representing the order in which
 * that subcell is visited by the Hilbert curve, plus 2 bits indicating the
 * new orientation of the Hilbert curve within that subcell. (Cell
 * orientations are represented as combination of kSwapMask and kInvertMask.)
 *
 * "lookup_ij" is an inverted table used for mapping in the opposite
 * direction.
 *
 * We also experimented with looking up 16 bits at a time (14 bits of position
 * plus 2 of orientation) but found that smaller lookup tables gave better
 * performance. (2KB fits easily in the primary cache.)
 * Values for these constants are *declared* in the *.h file. Even though
 * the declaration specifies a value for the constant, that declaration
 * is not a *definition* of storage for the value. Because the values are
 * supplied in the declaration, we don't need the values here. Failing to
 * define storage causes link errors for any code that tries to take the
 * address of one of these values.
 */
@ExperimentalUnsignedTypes
object LookupCellTables {

    const val kLookupBits = 4
    val lookupPos = UIntArray(1 shl 2 * kLookupBits + 2)
    val lookupIj = UIntArray(1 shl 2 * kLookupBits + 2)

    init {
        initLookupCell(0, 0, 0, 0, 0, 0)
        initLookupCell(0, 0, 0, S2Coords.kSwapMask, 0, S2Coords.kSwapMask)
        initLookupCell(0, 0, 0, S2Coords.kInvertMask, 0, S2Coords.kInvertMask)
        initLookupCell(0, 0, 0, S2Coords.kSwapMask or S2Coords.kInvertMask, 0, (S2Coords.kSwapMask or S2Coords.kInvertMask))
    }

    private fun initLookupCell(level: Int, i: Int, j: Int, origOrientation: Int, pos: Int, orientation: Int) {
        var currentLevel = level
        var currentI = i
        var currentJ = j
        var currentPos = pos
        if (currentLevel == kLookupBits) {
            val ij = (currentI shl kLookupBits) + currentJ
            lookupPos[(ij shl 2) + origOrientation] = ((currentPos shl 2) + orientation).toUInt()
            lookupIj[(currentPos shl 2) + origOrientation] = ((ij shl 2) + orientation).toUInt()
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

    /*
    fun getBits(n: LongArray, i: Int, j: Int, k: Int, bits: Int): Int {
        var bits = bits
        val mask = (1 shl kLookupBits) - 1
        bits += i shr k * kLookupBits and mask shl kLookupBits + 2
        bits += j shr k * kLookupBits and mask shl 2
        bits = lookupPos[bits]
        n[k shr 2] = n[k shr 2] or (bits.toLong() shr 2 shl (k and 3) * 2 * kLookupBits)
        bits = bits and (kSwapMask or kInvertMask)
        return bits
    }*/


}