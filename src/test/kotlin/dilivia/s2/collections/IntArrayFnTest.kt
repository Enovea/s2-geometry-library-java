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
package dilivia.s2.collections

import junit.framework.TestCase

class IntArrayFnTest : TestCase() {

    fun testIntArrayLowerBound() {
        val vector = intArrayOf(1,2,3,4,5,6)
        assertEquals(2, vector.lowerBound(1, 4, 3))
        assertEquals(4, vector.lowerBound(1, 4, 7))
        assertEquals(4, vector.lowerBound(1, 4, 8))
        assertEquals(2, vector.lowerBound(1, vector.size, 3))
        assertEquals(vector.size, vector.lowerBound(1, vector.size, 8))
        assertEquals(4, vector.lowerBound(4, vector.size, 3))
        assertEquals(vector.size, vector.lowerBound(vector.size + 1, vector.size + 3, 3))
    }

    fun testIntArrayLowerBoundEmpty() {
        val vector = intArrayOf()
        assertEquals(0, vector.lowerBound(1, 4, 3))
        assertEquals(0, vector.lowerBound(1, 4, 7))
        assertEquals(0, vector.lowerBound(1, 4, 8))
        assertEquals(0, vector.lowerBound(1, vector.size, 3))
        assertEquals(vector.size, vector.lowerBound(1, vector.size, 8))
        assertEquals(0, vector.lowerBound(4, vector.size, 3))
    }

    fun testIntArrayUpperBound() {
        val vector = intArrayOf(1,2,3,4,5,6)
        assertEquals(3, vector.upperBound(1, 4, 3))
        assertEquals(4, vector.upperBound(1, 4, 7))
        assertEquals(4, vector.upperBound(1, 4, 8))
        assertEquals(4, vector.upperBound(1, vector.size, 4))
        assertEquals(vector.size, vector.upperBound(1, vector.size, 8))
        assertEquals(4, vector.upperBound(4, vector.size, 3))
        assertEquals(vector.size, vector.upperBound(vector.size + 1, vector.size + 3, 3))
    }

    fun testIntArrayUpperBoundEmpty() {
        val vector = intArrayOf()
        assertEquals(0, vector.upperBound(1, 4, 3))
        assertEquals(0, vector.upperBound(1, 4, 7))
        assertEquals(0, vector.upperBound(1, 4, 8))
        assertEquals(0, vector.upperBound(1, vector.size, 3))
        assertEquals(vector.size, vector.upperBound(1, vector.size, 8))
        assertEquals(0, vector.upperBound(4, vector.size, 3))
    }
}
