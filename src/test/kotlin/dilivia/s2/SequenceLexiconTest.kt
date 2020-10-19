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

class SequenceLexiconTest : S2GeometryTestCase() {

    fun testint64() {
        val lex = SequenceLexicon<Long>()
        assertEquals(0, lex.add(listOf()))
        assertEquals(1, lex.add(listOf( 5L )))
        assertEquals(0, lex.add(listOf()))
        assertEquals(2, lex.add(listOf( 5L, 5L )))
        assertEquals(3, lex.add(listOf( 5L, 0L, -3L )))
        assertEquals(1, lex.add(listOf( 5L )))
        assertEquals(4, lex.add(listOf( 0x7fffffffffffffffL )))
        assertEquals(3, lex.add(listOf( 5, 0, -3 )))
        assertEquals(0, lex.add(listOf()))
        assertEquals(5, lex.size())
        expectSequence(listOf(), lex.sequence(0))
        expectSequence(listOf( 5L ), lex.sequence(1))
        expectSequence(listOf( 5L, 5L ), lex.sequence(2))
        expectSequence(listOf( 5L, 0L, -3L ), lex.sequence(3))
        expectSequence(listOf( 0x7fffffffffffffffL ), lex.sequence(4))
    }

    fun testClear() {
        val lex = SequenceLexicon<Long>()
        assertEquals(0, lex.add(listOf( 1L )))
        assertEquals(1, lex.add(listOf( 2L )))
        lex.clear()
        assertEquals(0, lex.add(listOf( 2L )))
        assertEquals(1, lex.add(listOf( 1L )))
    }

    companion object {
        fun <T> expectSequence(expected: List<T>, actual: SequenceLexicon<T>.Sequence) {
            assertEquals(expected.size, actual.size())
            assertEquals(expected, actual.values())
        }
    }

}

