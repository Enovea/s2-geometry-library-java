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

import dilivia.s2.S2Coords.kIJtoPos
import dilivia.s2.S2Coords.kInvertMask
import dilivia.s2.S2Coords.kPosToIJ
import dilivia.s2.S2Coords.kSwapMask
import dilivia.s2.S2Random.randomCellId
import dilivia.s2.S2Random.randomInt
import kotlin.math.abs

@ExperimentalUnsignedTypes
class S2CoordsTest : S2GeometryTestCase() {

    fun swapAxes(ij: Int): Int = ((ij shr 1) and 1) + ((ij and 1) shl 1)

    fun invertBits(ij: Int): Int = ij xor 3

    fun testTraversalOrder() {
        for (r in 0..3) {
            for (i in 0..3) {
                // Check consistency with respect to swapping axes.
                assertEquals(kIJtoPos[r][i], kIJtoPos[r xor kSwapMask][swapAxes(i)])
                assertEquals(kPosToIJ[r][i], swapAxes(kPosToIJ[r xor kSwapMask][i]))

                // Check consistency with respect to reversing axis directions.
                assertEquals(kIJtoPos[r][i], kIJtoPos[r xor kInvertMask][invertBits(i)])
                assertEquals(kPosToIJ[r][i], invertBits(kPosToIJ[r xor kInvertMask][i]))

                // Check that the two tables are inverses of each other.
                assertEquals(kIJtoPos[r][kPosToIJ[r][i]], i)
                assertEquals(kPosToIJ[r][kIJtoPos[r][i]], i)
            }
        }
    }

    fun testST_UV_Conversions() {
        // Check boundary conditions.
        var s = 0.0
        while (s <= 1) {
            val u = S2Coords.stToUv(s)
            assertEquals(u, 2 * s - 1)
            s += 0.5
        }
        for (u in -1..1) {
            s = S2Coords.uvToSt(u.toDouble())
            assertEquals(s, 0.5 * (u + 1))
        }
        // Check that UVtoST and STtoUV are inverses.
        var x = 0.0
        while (x <= 1.0) {
            assertEquals(S2Coords.uvToSt(S2Coords.stToUv(x)), x, 1e-15)
            assertEquals(S2Coords.stToUv(S2Coords.uvToSt(2 * x - 1)), 2 * x - 1, 1e-15)
            x += 0.0001
        }
    }

    fun testFaceUVtoXYZ() {
        // Check that each face appears exactly once.
        var sum = S2Point()
        for (face in 0..5) {
            val center = S2Coords.faceUVtoXYZ(face, 0.0, 0.0)
            assertEquals(S2Coords.getNorm(face), center)
            assertEquals(abs(center[center.largestAbsComponent()]), 1.0)
            sum += center.abs()
        }
        assertEquals(sum, S2Point(2, 2, 2))

        // Check that each face has a right-handed coordinate system.
        for (face in 0..5) {
            assertEquals(S2Coords.getUAxis(face)
                    .crossProd(S2Coords.getVAxis(face))
                    .dotProd(S2Coords.faceUVtoXYZ(face, 0.0, 0.0)), 1.0)
        }

        // Check that the Hilbert curves on each face combine to form a
        // continuous curve over the entire cube.
        for (face in 0..5) {
            // The Hilbert curve on each face starts at (-1,-1) and terminates
            // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
            val sign = if (face and kSwapMask != 0) -1.0 else 1.0
            assertEquals(S2Coords.faceUVtoXYZ(face, sign, -sign), S2Coords.faceUVtoXYZ((face + 1) % 6, -1.0, -1.0))
        }
    }

    fun testFaceXYZtoUVW() {
        for (face in 0..5) {
            assertEquals(S2Point(0, 0, 0), S2Coords.faceXYZtoUVW(face, S2Point(0, 0, 0)))
            assertEquals(S2Point(1, 0, 0), S2Coords.faceXYZtoUVW(face, S2Coords.getUAxis(face)))
            assertEquals(S2Point(-1, 0, 0), S2Coords.faceXYZtoUVW(face, -S2Coords.getUAxis(face)))
            assertEquals(S2Point(0, 1, 0), S2Coords.faceXYZtoUVW(face, S2Coords.getVAxis(face)))
            assertEquals(S2Point(0, -1, 0), S2Coords.faceXYZtoUVW(face, -S2Coords.getVAxis(face)))
            assertEquals(S2Point(0, 0, 1), S2Coords.faceXYZtoUVW(face, S2Coords.getNorm(face)))
            assertEquals(S2Point(0, 0, -1), S2Coords.faceXYZtoUVW(face, -S2Coords.getNorm(face)))
        }
    }

    fun testXYZToFaceSiTi() {
        // Check the conversion of random cells to center points and back.
        for (level in 0..S2CellId.kMaxLevel) {
            for (i in 0 until 1000) {
                val id = randomCellId(level)

                val (actual_level, faceSiTi) = S2Coords.xyzToFaceSiTi(id.toPoint())
                assertEquals(level, actual_level)
                val actual_id = S2CellId.fromFaceIJ(faceSiTi.face, (faceSiTi.si / 2U).toInt(), (faceSiTi.ti / 2U).toInt()).parent(level)
                assertEquals(id, actual_id)

                // Now test a point near the cell center but not equal to it.
                val p_moved = id.toPoint() + S2Point(1e-13, 1e-13, 1e-13)
                val (actual_level_moved, faceSiTi_moved) = S2Coords.xyzToFaceSiTi(p_moved)
                assertEquals(-1, actual_level_moved)
                assertEquals(faceSiTi, faceSiTi_moved)

                // Finally, test some random (si,ti) values that may be at different
                // levels, or not at a valid level at all (for example, si == 0).
                val face_random = randomInt(S2CellId.kNumFaces)
                var si_random = 0U
                var ti_random = 0U
                val mask = -1 shl (S2CellId.kMaxLevel - level)
                do {
                    si_random = (randomInt() and mask).toUInt()
                    ti_random = (randomInt() and mask).toUInt()
                } while (si_random > S2Coords.kMaxSiTi || ti_random > S2Coords.kMaxSiTi)
                val p_random = S2Coords.faceSiTitoXYZ(face_random, si_random, ti_random)
                val (actual_level_random, faceSiTi_random) = S2Coords.xyzToFaceSiTi(p_random)
                if (faceSiTi_random.face != face_random) {
                    // The chosen point is on the edge of a top-level face cell.
                    assertEquals(-1, actual_level_random)
                    assertTrue(faceSiTi_random.si == 0U || faceSiTi_random.si == S2Coords.kMaxSiTi
                            || faceSiTi_random.ti == 0U || faceSiTi_random.ti == S2Coords.kMaxSiTi)
                } else {
                    assertEquals(si_random, faceSiTi_random.si)
                    assertEquals(ti_random, faceSiTi_random.ti)
                    if (actual_level_random >= 0) {
                        assertEquals(p_random, S2CellId.fromFaceIJ(faceSiTi_random.face, (faceSiTi_random.si / 2U).toInt(), (faceSiTi_random.ti / 2U).toInt()).parent(actual_level_random).toPoint())
                    }
                }
            }
        }
    }

    fun testUVNorms() {
        // Check that GetUNorm and GetVNorm compute right-handed normals for
        // an edge in the increasing U or V direction.
        for (face in 0..5) {
            var x = -1.0
            while (x <= 1.0) {
                assertEquals(S2Coords.faceUVtoXYZ(face, x, -1.0).crossProd(S2Coords.faceUVtoXYZ(face, x, 1.0)).angle(S2Coords.getUNorm(face, x)), 0.0)
                assertEquals(S2Coords.faceUVtoXYZ(face, -1.0, x).crossProd(S2Coords.faceUVtoXYZ(face, 1.0, x)).angle(S2Coords.getVNorm(face, x)), 0.0)
                x += 1.0 / 1024.0
            }
        }
    }

    fun testUVWAxis() {
        for (face in 0..5) {
            // Check that axes are consistent with FaceUVtoXYZ.
            assertEquals(S2Coords.faceUVtoXYZ(face, 1.0, 0.0) - S2Coords.faceUVtoXYZ(face, 0.0, 0.0), S2Coords.getUAxis(face))
            assertEquals(S2Coords.faceUVtoXYZ(face, 0.0, 1.0) - S2Coords.faceUVtoXYZ(face, 0.0, 0.0), S2Coords.getVAxis(face))
            assertEquals(S2Coords.faceUVtoXYZ(face, 0.0, 0.0), S2Coords.getNorm(face))

            // Check that every face coordinate frame is right-handed.
            assertEquals(1.0, S2Coords.getUAxis(face).crossProd(S2Coords.getVAxis(face)).dotProd(S2Coords.getNorm(face)))

            // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
            assertEquals(S2Coords.getUAxis(face), S2Coords.getUVWAxis(face, 0))
            assertEquals(S2Coords.getVAxis(face), S2Coords.getUVWAxis(face, 1))
            assertEquals(S2Coords.getNorm(face), S2Coords.getUVWAxis(face, 2))
        }
    }

    fun testUVWFace() {
        // Check that GetUVWFace is consistent with GetUVWAxis.
        for (face in 0..5) {
            for (axis in 0..2) {
                assertEquals(S2Coords.getFace(-S2Coords.getUVWAxis(face, axis)), S2Coords.getUVWFace(face, axis, 0))
                assertEquals(S2Coords.getFace(S2Coords.getUVWAxis(face, axis)), S2Coords.getUVWFace(face, axis, 1))
            }
        }
    }

}