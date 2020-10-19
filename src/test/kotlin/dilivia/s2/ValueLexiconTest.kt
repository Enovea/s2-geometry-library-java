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

class ValueLexiconTest : S2GeometryTestCase() {

  fun testDuplicateValues() {
    val lex = ValueLexicon<Long>()
    assertEquals(0, lex.add(5L))
    assertEquals(1, lex.add(0L))
    assertEquals(1, lex.add(0L))
    assertEquals(2, lex.add(-3L))
    assertEquals(0, lex.add(5L))
    assertEquals(1, lex.add(0L))
    assertEquals(3, lex.add(0x7fffffffffffffffL))
    assertEquals(4, lex.add(-0x80L shl (14*4)))
    assertEquals(3, lex.add(0x7fffffffffffffffL))
    assertEquals(4, lex.add(-0x80L shl (14*4)))
    assertEquals(5, lex.size())
    assertEquals(5L, lex.value(0))
    assertEquals(0L, lex.value(1))
    assertEquals(-3L, lex.value(2))
    assertEquals(0x7fffffffffffffff, lex.value(3))
    assertEquals(-0x80L shl (14*4), lex.value(4))
  }

  fun testClear() {
    val lex = ValueLexicon<Long>()
    assertEquals(0, lex.add(1))
    assertEquals(1, lex.add(2))
    assertEquals(0, lex.add(1))
    lex.clear()
    assertEquals(0, lex.add(2))
    assertEquals(1, lex.add(1))
    assertEquals(0, lex.add(2))
  }

  fun testFloatEquality() {
    val lex = ValueLexicon<S2Point>()
    val a = S2Point (1.0, 0.0, 0.0)
    val b = S2Point (1.0, -0.0, 0.0)
    val c = S2Point(1.0, 0.0, -0.0)
    assertEquals(0, lex.add(a))
    assertEquals(0, lex.add(b))
    assertEquals(0, lex.add(c))
    assertEquals(1, lex.size())
  }

}
