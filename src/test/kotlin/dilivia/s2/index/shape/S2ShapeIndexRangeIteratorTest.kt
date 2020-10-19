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

import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser

class S2ShapeIndexRangeIteratorTest : S2GeometryTestCase() {

    fun testNext() {
        // Create an index with one point each on S2CellId faces 0, 1, and 2.
        val index: S2ShapeIndex = S2TextParser.makeIndex("0:0 | 0:90 | 90:0 # #")
        val it = S2ShapeIndexRangeIterator(index)
        assertEquals(0, it.id().face())
        it.next()
        assertEquals(1, it.id().face())
        it.next()
        assertEquals(2, it.id().face())
        it.next()
        assertEquals(S2CellId.sentinel(), it.id())
        assertTrue(it.done())
    }

    fun testEmptyIndex() {
        val empty: S2ShapeIndex = S2TextParser.makeIndex("# #")
        val non_empty: S2ShapeIndex = S2TextParser.makeIndex("0:0 # #")
        val empty_it = S2ShapeIndexRangeIterator(empty)
        val non_empty_it = S2ShapeIndexRangeIterator(non_empty)
        assertFalse(non_empty_it.done())
        assertTrue(empty_it.done())

        empty_it.seekTo(non_empty_it)
        assertTrue(empty_it.done())

        empty_it.seekBeyond(non_empty_it)
        assertTrue(empty_it.done())

        empty_it.seekTo(empty_it)
        assertTrue(empty_it.done())

        empty_it.seekBeyond(empty_it)
        assertTrue(empty_it.done())
    }

} 
