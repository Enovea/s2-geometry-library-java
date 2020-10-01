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
package dilivia.s2.shape

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2Random

class S2EdgeVectorShapeTest : S2GeometryTestCase() {

    fun testEmpty() {
        val shape = S2EdgeVectorShape()
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(1, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testEdgeAccess() {
        val shape = S2EdgeVectorShape()
        S2Random.reset(0)
        val kNumEdges = 100
        repeat(kNumEdges) {
            val a = S2Random.randomPoint()  // Control the evaluation order
            shape.add(a, S2Random.randomPoint())
        }
        assertEquals(kNumEdges, shape.numEdges)
        assertEquals(kNumEdges, shape.numChains)
        assertEquals(1, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        S2Random.reset(0)
        repeat(kNumEdges) { i ->
            assertEquals(i, shape.chain(i).start)
            assertEquals(1, shape.chain(i).length)
            val edge = shape.edge(i)
            assertEquals(S2Random.randomPoint(), edge.v0)
            assertEquals(S2Random.randomPoint(), edge.v1)
        }
    }

    fun testSingletonConstructor() {
        val a = S2Point(1, 0, 0)
        val b = S2Point(0, 1, 0)
        val shape = S2EdgeVectorShape(a, b)
        assertEquals(1, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        val edge = shape.edge(0)
        assertEquals(a, edge.v0)
        assertEquals(b, edge.v1)
    }
}