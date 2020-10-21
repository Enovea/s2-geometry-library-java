/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
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
package dilivia.s2.shape

import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.shape.S2Shape.Companion.containsBruteForce

class GetReferencePointTest : S2GeometryTestCase() {

    fun testEmptyPolygon() {
        val shape = S2LaxPolygonShape(S2Polygon())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testFullPolygon() {
        val shape = S2LaxPolygonShape(S2Polygon(S2TextParser.makeLoop("full")))
        assertTrue(shape.getReferencePoint().contained)
    }

    fun testDegenerateLoops() {
        val loops = listOf(
                S2TextParser.parsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
                S2TextParser.parsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
                S2TextParser.parsePoints("5:5, 6:6")
        )
        val shape = S2LaxPolygonShape(loops)
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testInvertedLoops() {
        val loops = listOf(
                S2TextParser.parsePoints("1:2, 1:1, 2:2"),
                S2TextParser.parsePoints("3:4, 3:3, 4:4")
        )
        val shape = S2LaxPolygonShape(loops)
        assertTrue(containsBruteForce(shape, S2Point.origin()))
    }

    fun testPartiallyDegenerateLoops() {
        repeat(100) {
            // First we construct a long convoluted edge chain that follows the
            // S2CellId Hilbert curve.  At some random point along the curve, we
            // insert a small triangular loop.
            val loops = mutableListOf<Loop>()
            val loop = mutableListOf<S2Point>()
            loops.add(loop)
            val num_vertices = 100L
            val start = S2Random.randomCellId(S2CellId.kMaxLevel - 1)
            val end = start.advanceWrap(num_vertices)
            val loop_cellid = start.advanceWrap(S2Random.randomLong(num_vertices - 2) + 1)
            val triangle = mutableListOf<S2Point>()
            var cellid = start
            while (cellid != end) {
                if (cellid == loop_cellid) {
                    // Insert a small triangular loop.  We save the loop so that we can
                    // test whether it contains the origin later.
                    triangle.add(cellid.child(0).toPoint())
                    triangle.add(cellid.child(1).toPoint())
                    triangle.add(cellid.child(2).toPoint())
                    loop.addAll(triangle)
                    loop.add(cellid.child(0).toPoint())
                } else {
                    loop.add(cellid.toPoint())
                }
                cellid = cellid.nextWrap()
            }
            // Now we retrace our steps, except that we skip the three edges that form
            // the triangular loop above.
            cellid = end
            while (cellid != start) {
                if (cellid == loop_cellid) {
                    loop.add(cellid.child(0).toPoint())
                } else {
                    loop.add(cellid.toPoint())
                }
                cellid = cellid.prevWrap()
            }
            val shape = S2LaxPolygonShape(loops)
            val triangle_loop = S2Loop(triangle)
            val ref = shape.getReferencePoint()
            assertEquals(triangle_loop.contains(ref.point), ref.contained)
        }
    }

}  
