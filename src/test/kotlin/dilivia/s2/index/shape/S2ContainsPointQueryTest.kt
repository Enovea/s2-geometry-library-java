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
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.S2TextParser
import dilivia.s2.index.shape.S2ContainsPointQuery.Companion.makeS2ContainsPointQuery
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Loop
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId

typealias EdgeIdVector = MutableList<ShapeEdgeId>

class S2ContainsPointQueryTest : S2GeometryTestCase() {

    fun testS2ContainsPointQueryVertexModelOpen() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.OPEN)
        val q = makeS2ContainsPointQuery(index, options)
        assertFalse(q.contains(S2TextParser.makePoint("0:0")))
        assertFalse(q.contains(S2TextParser.makePoint("-1:1")))
        assertFalse(q.contains(S2TextParser.makePoint("1:1")))
        assertFalse(q.contains(S2TextParser.makePoint("0:2")))
        assertFalse(q.contains(S2TextParser.makePoint("0:3")))
        assertFalse(q.contains(S2TextParser.makePoint("0:5")))
        assertFalse(q.contains(S2TextParser.makePoint("0:7")))
        assertFalse(q.contains(S2TextParser.makePoint("2:6")))
        assertTrue(q.contains(S2TextParser.makePoint("1:6")))
        assertFalse(q.contains(S2TextParser.makePoint("10:10")))

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertFalse(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6")))
        assertFalse(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5")))
        assertFalse(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7")))
    }

    fun testS2ContainsPointQueryVertexModelSemiOpen() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.SEMI_OPEN)
        val q = makeS2ContainsPointQuery(index, options)
        assertFalse(q.contains(S2TextParser.makePoint("0:0")))
        assertFalse(q.contains(S2TextParser.makePoint("-1:1")))
        assertFalse(q.contains(S2TextParser.makePoint("1:1")))
        assertFalse(q.contains(S2TextParser.makePoint("0:2")))
        assertFalse(q.contains(S2TextParser.makePoint("0:5")))
        assertTrue(q.contains(S2TextParser.makePoint("0:7")))  // Contained vertex.
        assertFalse(q.contains(S2TextParser.makePoint("2:6")))
        assertTrue(q.contains(S2TextParser.makePoint("1:6")))
        assertFalse(q.contains(S2TextParser.makePoint("10:10")))

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertFalse(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6")))
        assertFalse(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7")))
    }

    fun testS2ContainsPointQueryVertexModelClosed() {
        val index = S2TextParser.makeIndex("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6")
        val options = S2ContainsPointQueryOptions(S2VertexModel.CLOSED)
        val q = makeS2ContainsPointQuery(index, options)
        assertTrue(q.contains(S2TextParser.makePoint("0:0")))
        assertTrue(q.contains(S2TextParser.makePoint("-1:1")))
        assertTrue(q.contains(S2TextParser.makePoint("1:1")))
        assertFalse(q.contains(S2TextParser.makePoint("0:2")))
        assertTrue(q.contains(S2TextParser.makePoint("0:5")))
        assertTrue(q.contains(S2TextParser.makePoint("0:7")))
        assertTrue(q.contains(S2TextParser.makePoint("2:6")))
        assertTrue(q.contains(S2TextParser.makePoint("1:6")))
        assertFalse(q.contains(S2TextParser.makePoint("10:10")))

        // Test the last few cases using the Init() method instead.
        val q2 = S2ContainsPointQuery<MutableS2ShapeIndex>()
        q2.init(index, options)
        assertFalse(q2.shapeContains(index.shape(1)!!, S2TextParser.makePoint("1:6")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("1:6")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:5")))
        assertTrue(q2.shapeContains(index.shape(2)!!, S2TextParser.makePoint("0:7")))
    }

    fun testS2ContainsPointQueryGetContainingShapes() {
        // Also tests shapeContains().
        val kNumVerticesPerLoop = 10
        val kMaxLoopRadius = kmToAngle(10.0)
        val centerCap = S2Cap.fromCenterAngle(S2Random.randomPoint(), kMaxLoopRadius)
        val index = MutableS2ShapeIndex()
        repeat(100) {
            val loop = S2Loop.makeRegularLoop(
                    S2Random.samplePoint(centerCap),
                    kMaxLoopRadius * S2Random.randomDouble(), kNumVerticesPerLoop)
            index.add(S2Loop.Shape(loop = loop))
        }
        val query = makeS2ContainsPointQuery(index)
        repeat(100) {
            val p = S2Random.samplePoint(centerCap)
            val expected = mutableListOf<S2Shape>()
            for (shape in index) {
                val loop = (shape as S2Loop.Shape).loop
                if (loop.contains(p)) {
                    assertTrue(query.shapeContains(shape, p))
                    expected.add(shape)
                } else {
                    assertFalse(query.shapeContains(shape, p))
                }
            }
            val actual = query.getContainingShapes(p)
            assertEquals(expected, actual)
        }
    }

    fun expectIncidentEdgeIds(expected: EdgeIdVector, index: MutableS2ShapeIndex, p: S2Point) {
        val actual: EdgeIdVector = mutableListOf()
        val q = makeS2ContainsPointQuery(index)
        assertTrue(q.visitIncidentEdges(p) { e ->
          actual.add(e.id)
          true
        })
        assertEquals(expected, actual)
    }

    fun testS2ContainsPointQueryVisitIncidentEdges() {
        val index = S2TextParser.makeIndex("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2")
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(0, 0)), index, S2TextParser.makePoint("0:0"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(0, 1), ShapeEdgeId(1, 0)), index, S2TextParser.makePoint("1:1"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(1, 0), ShapeEdgeId(2, 0), ShapeEdgeId(2, 2)), index, S2TextParser.makePoint("1:2"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(2, 0), ShapeEdgeId(2, 1)), index, S2TextParser.makePoint("1:3"))
        expectIncidentEdgeIds(mutableListOf(ShapeEdgeId(2, 1), ShapeEdgeId(2, 2)), index, S2TextParser.makePoint("2:2"))
    }

}
