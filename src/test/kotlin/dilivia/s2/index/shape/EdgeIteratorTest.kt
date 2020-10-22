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
import dilivia.s2.shape.Edge

class EdgeIteratorTest : S2GeometryTestCase() {

    companion object {

        // Returns the full list of edges in g.
        // The edges are collected from points, lines, and polygons in that order.
        fun getEdges(index: S2ShapeIndex): List<Edge> {
            val result = mutableListOf<Edge>()
            for (shape in index) {
                if (shape == null) continue
                for (j in 0 until shape.numEdges) {
                    result.add(shape.edge(j))
                }
            }
            return result
        }

        // Verifies that the edges produced by an EdgeIterator matches GetEdges.
        fun verify(index: S2ShapeIndex) {
            val expected = getEdges(index);

            var i = 0;
            val it = EdgeIterator(index)
            while (!it.done()) {
                assertTrue(i < expected.size)
                assertEquals(expected[i], it.edge())
                it.next()
                ++i
            }
        }

    }

    fun testEmpty() {
        val index = S2TextParser.makeIndex("##")
        verify(index)
    }

    fun testPoints() {
        val index = S2TextParser.makeIndex("0:0|1:1##")
        verify(index)
    }

    fun testLines() {
        val index = S2TextParser.makeIndex("#0:0,10:10|5:5,5:10|1:2,2:1#")
        verify(index)
    }

    fun testPolygons() {
        val index = S2TextParser.makeIndex("##10:10,10:0,0:0|-10:-10,-10:0,0:0,0:-10")
        verify(index)
    }

    fun testCollection() {
        val index = S2TextParser.makeIndex(
                "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0")
        verify(index)
    }

    fun testRemove() {
        val index = S2TextParser.makeIndex(
                "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0")
        index.remove(0)
        verify(index)
    }

    fun testAssignmentAndEquality() {
        val index1 = S2TextParser.makeIndex(
                "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

        val index2 = S2TextParser.makeIndex(
                "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

        var it1 = EdgeIterator(index1)
        val it2 = EdgeIterator(index2)

        // Different indices.
        assertTrue(it1 != it2);

        it1 = it2.clone()
        assertEquals(it1, it2)

        it1.next();
        assertTrue(it1 != it2);

        it2.next();
        assertEquals(it1, it2)
    }

} 
