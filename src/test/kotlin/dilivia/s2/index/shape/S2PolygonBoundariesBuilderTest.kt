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
package dilivia.s2.index.shape

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser
import dilivia.s2.shape.S2LaxLoopShape
import dilivia.s2.shape.S2Shape


class TestLaxLoop(vertex_str: String) : S2LaxLoopShape() {
    init {
        val vertices = S2TextParser.parsePoints(vertex_str)
        init(vertices)
    }
}

class S2PolygonBoundariesBuilderTest : S2GeometryTestCase() {

    fun testBuildPolygonBoundariesNoComponents() {
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>()
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(0, faces.size)
    }

    fun testBuildPolygonBoundariesOneLoop() {
        val a0 = TestLaxLoop("0:0, 1:0, 0:1")  // Outer face
        val a1 = TestLaxLoop("0:0, 0:1, 1:0")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(listOf(a0, a1))
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(2, faces.size);
    }

    fun testBuildPolygonBoundariesTwoLoopsSameComponent() {
        val a0 = TestLaxLoop("0:0, 1:0, 0:1")  // Outer face
        val a1 = TestLaxLoop("0:0, 0:1, 1:0")
        val a2 = TestLaxLoop("1:0, 0:1, 1:1")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(listOf(a0, a1, a2))
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(3, faces.size)
    }

    fun testBuildPolygonBoundariesTwoNestedLoops() {
        val a0 = TestLaxLoop("0:0, 3:0, 0:3")  // Outer face
        val a1 = TestLaxLoop("0:0, 0:3, 3:0")
        val b0 = TestLaxLoop("1:1, 2:0, 0:2")  // Outer face
        val b1 = TestLaxLoop("1:1, 0:2, 2:0")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(
                listOf(a0, a1),
                listOf(b0, b1)
        )
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(3, faces.size);
        assertEquals(mutableListOf(b0, a1), faces[0])
    }

    fun testBuildPolygonBoundariesTwoLoopsDifferentComponents() {
        val a0 = TestLaxLoop("0:0, 1:0, 0:1")  // Outer face
        val a1 = TestLaxLoop("0:0, 0:1, 1:0")
        val b0 = TestLaxLoop("0:2, 1:2, 0:3")  // Outer face
        val b1 = TestLaxLoop("0:2, 0:3, 1:2")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(
                listOf(a0, a1),
                listOf(b0, b1)
        )
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(3, faces.size)
        assertEquals((mutableListOf<S2Shape>(a0, b0)), faces[2])
    }

    fun testBuildPolygonBoundariesOneDegenerateLoop() {
        val a0 = TestLaxLoop("0:0, 1:0, 0:0")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(listOf(a0))
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(1, faces.size)
    }

    fun testBuildPolygonBoundariesTwoDegenerateLoops() {
        val a0 = TestLaxLoop("0:0, 1:0, 0:0")
        val b0 = TestLaxLoop("2:0, 3:0, 2:0")
        val faces = mutableListOf<MutableList<S2Shape>>()
        val components = listOf<List<S2Shape>>(
                listOf(a0),
                listOf(b0)
        )
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(1, faces.size)
        assertEquals(2, faces[0].size)
    }

    fun testBuildPolygonBoundariesComplexTest1() {
        // Loops at index 0 are the outer (clockwise) loops.
        // Component "a" consists of 4 adjacent squares forming a larger square.
        val a0 = TestLaxLoop("0:0, 25:0, 50:0, 50:25, 50:50, 25:50, 0:50, 0:50")
        val a1 = TestLaxLoop("0:0, 0:25, 25:25, 25:0")
        val a2 = TestLaxLoop("0:25, 0:50, 25:50, 25:25")
        val a3 = TestLaxLoop("25:0, 25:25, 50:25, 50:0")
        val a4 = TestLaxLoop("25:25, 25:50, 50:50, 50:25")
        // Component "b" consists of a degenerate loop to the left of "a".
        val b0 = TestLaxLoop("0:-10, 10:-10")
        // Components "a1_a", "a1_b", and "a1_c" are located within "a1".
        val a1_a0 = TestLaxLoop("5:5, 20:5, 20:10, 5:10")
        val a1_a1 = TestLaxLoop("5:5, 5:10, 10:10, 10:5")
        val a1_a2 = TestLaxLoop("10:5, 10:10, 15:10, 15:5")
        val a1_a3 = TestLaxLoop("15:5, 15:10, 20:10, 20:5")
        val a1_b0 = TestLaxLoop("5:15, 20:15, 20:20, 5:20")
        val a1_b1 = TestLaxLoop("5:15, 5:20, 20:20, 20:15")
        val a1_c0 = TestLaxLoop("2:5, 2:10, 2:5")
        // Two components located inside "a1_a2" and "a1_a3".
        val a1_a2_a0 = TestLaxLoop("11:6, 14:6, 14:9, 11:9")
        val a1_a2_a1 = TestLaxLoop("11:6, 11:9, 14:9, 14:6")
        val a1_a3_a0 = TestLaxLoop("16:6, 19:9, 16:6")
        // Five component located inside "a3" and "a4".
        val a3_a0 = TestLaxLoop("30:5, 45:5, 45:20, 30:20")
        val a3_a1 = TestLaxLoop("30:5, 30:20, 45:20, 45:5")
        val a4_a0 = TestLaxLoop("30:30, 40:30, 30:30")
        val a4_b0 = TestLaxLoop("30:35, 40:35, 30:35")
        val a4_c0 = TestLaxLoop("30:40, 40:40, 30:40")
        val a4_d0 = TestLaxLoop("30:45, 40:45, 30:45")
        val components = listOf(
                listOf(a0, a1, a2, a3, a4),
                listOf(b0),
                listOf(a1_a0, a1_a1, a1_a2, a1_a3),
                listOf(a1_b0, a1_b1),
                listOf(a1_c0),
                listOf(a1_a2_a0, a1_a2_a1),
                listOf(a1_a3_a0),
                listOf(a3_a0, a3_a1),
                listOf(a4_a0),
                listOf(a4_b0),
                listOf(a4_c0),
                listOf(a4_d0)
        )
        val expected_faces = listOf(
                listOf(a0, b0),
                listOf(a1, a1_a0, a1_b0, a1_c0),
                listOf(a1_a1),
                listOf(a1_a2, a1_a2_a0),
                listOf(a1_a2_a1),
                listOf(a1_a3, a1_a3_a0),
                listOf(a1_b1),
                listOf(a2),
                listOf(a3, a3_a0),
                listOf(a3_a1),
                listOf(a4, a4_a0, a4_b0, a4_c0, a4_d0)
        )
        val faces = mutableListOf<MutableList<S2Shape>>()
        S2PolygonBoundariesBuilder.buildPolygonBoundaries(components, faces)
        assertEquals(expected_faces.size, faces.size)
        assertEquals(expected_faces, faces)
    }

}
