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
import dilivia.s2.S2TextParser

class S2LaxPolylineShapeTest : S2GeometryTestCase() {

    fun testNoVertices() {
        val shape = S2LaxPolylineShape()
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(1, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testOneVertex() {
        val vertices = listOf(S2Point(1, 0, 0))
        val shape = S2LaxPolylineShape(vertices)
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(1, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
    }

    fun testEdgeAccess() {
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1")
        val shape = S2LaxPolylineShape(vertices)
        assertEquals(2, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(2, shape.chain(0).length)
        assertEquals(1, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        val edge0 = shape.edge(0)
        assertEquals(vertices[0], edge0.v0)
        assertEquals(vertices[1], edge0.v1)
        val edge1 = shape.edge(1)
        assertEquals(vertices[1], edge1.v0)
        assertEquals(vertices[2], edge1.v1)
    }

}