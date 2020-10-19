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
import dilivia.s2.S2TextParser
import dilivia.s2.region.S2Loop

class S2LaxLoopShapeTest : S2GeometryTestCase() {

    fun testS2LaxLoopShapeEmptyLoop() {
        // Test S2Loop constructor.
        val shape = S2LaxLoopShape()
        shape.init(S2Loop(S2Loop.kEmpty))
        assertEquals(0, shape.numVertices())
        assertEquals(0, shape.numEdges)
        assertEquals(0, shape.numChains)
        assertEquals(2, shape.dimension)
        assertTrue(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testS2LaxLoopShapeNonEmptyLoop() {
        // Test vector<S2Point> constructor.
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val shape = S2LaxLoopShape(vertices)
        assertEquals(vertices.size, shape.numVertices())
        assertEquals(vertices.size, shape.numEdges)
        assertEquals(1, shape.numChains)
        assertEquals(0, shape.chain(0).start)
        assertEquals(vertices.size, shape.chain(0).length)
        for (i in vertices.indices) {
            assertEquals(vertices[i], shape.vertex(i))
            val edge = shape . edge (i)
            assertEquals(vertices[i], edge.v0)
            assertEquals(vertices[(i + 1) % vertices.size], edge.v1)
        }
        assertEquals(2, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }

    fun testS2LaxClosedPolylineShapeNoInterior() {
        val vertices = S2TextParser.parsePoints("0:0, 0:1, 1:1, 1:0")
        val shape = S2LaxClosedPolylineShape(vertices)
        assertEquals(1, shape.dimension)
        assertFalse(shape.isEmpty())
        assertFalse(shape.isFull())
        assertFalse(shape.getReferencePoint().contained)
    }
/*
  fun testS2VertexIdLaxLoopShapeEmptyLoop()
  {
    S2VertexIdLaxLoopShape shape (vector<int32>(), nullptr)
    assertEquals(0, shape.numEdges)
    assertEquals(0, shape.numVertices())
    assertEquals(0, shape.numChains)
    assertEquals(2, shape.dimension)
    assertTrue(shape.isEmpty())
    assertFalse(shape.isFull())
    assertFalse(shape.getReferencePoint().contained)
  }

  fun testS2VertexIdLaxLoopShapeInvertedLoop()
  {
    vector<S2Point> vertex_array =
    s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0")
    vector<int32> vertex_ids { 0, 3, 2, 1 };  // Inverted.
    S2VertexIdLaxLoopShape shape (vertex_ids, &vertex_array[0])
    assertEquals(4, shape.numEdges)
    assertEquals(4, shape.numVertices())
    assertEquals(1, shape.numChains)
    assertEquals(0, shape.chain(0).start)
    assertEquals(4, shape.chain(0).length)
    assertEquals(& vertex_array [0], &shape.vertex(0))
    assertEquals(& vertex_array [3], &shape.vertex(1))
    assertEquals(& vertex_array [2], &shape.vertex(2))
    assertEquals(& vertex_array [1], &shape.vertex(3))
    assertEquals(2, shape.dimension)
    assertFalse(shape.isEmpty())
    assertFalse(shape.isFull())
    assertTrue(s2shapeutil::ContainsBruteForce(shape, S2::Origin()))
  }
*/
}
