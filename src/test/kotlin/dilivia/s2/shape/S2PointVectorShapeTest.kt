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

class S2PointVectorShapeTest : S2GeometryTestCase() {

    fun testEmpty() {
        val shape = S2PointVectorShape()
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(0, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testConstructionAndAccess() {
        val points = mutableListOf<S2Point>()
        S2Random.reset(0);
        val kNumPoints = 100;
        repeat(kNumPoints) {
            points.add(S2Random.randomPoint())
        }
        val shape = S2PointVectorShape(points = points)

        assertEquals(kNumPoints, shape.numEdges)
        assertEquals(kNumPoints, shape.numChains)
        assertEquals(0, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        S2Random.reset(0)
        repeat(kNumPoints) { i ->
            assertEquals(i, shape.chain(i).start)
            assertEquals(1, shape.chain(i).length)
            val edge = shape.edge(i)
            val pt = S2Random.randomPoint()
            assertEquals(pt, edge.v0)
            assertEquals(pt, edge.v1)
        }
    }

}
