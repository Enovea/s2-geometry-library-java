package dilivia.s2.collections

import dilivia.s2.collections.lowerBound
import junit.framework.TestCase

class ListFnTest : TestCase() {


    fun testListLowerBound() {
        val vector = listOf(1,2,3,4,5,6)
        assertEquals(2, vector.lowerBound(1, 4, 3))
        assertEquals(4, vector.lowerBound(1, 4, 7))
        assertEquals(4, vector.lowerBound(1, 4, 8))
        assertEquals(2, vector.lowerBound(1, vector.size, 3))
        assertEquals(vector.size, vector.lowerBound(1, vector.size, 8))
        assertEquals(4, vector.lowerBound(4, vector.size, 3))
        assertEquals(vector.size, vector.lowerBound(vector.size + 1, vector.size + 3, 3))
    }

    fun testListLowerBoundEmpty() {
        val vector = listOf<Int>()
        assertEquals(0, vector.lowerBound(1, 4, 3))
        assertEquals(0, vector.lowerBound(1, 4, 7))
        assertEquals(0, vector.lowerBound(1, 4, 8))
        assertEquals(0, vector.lowerBound(1, vector.size, 3))
        assertEquals(0, vector.lowerBound(1, vector.size, 8))
        assertEquals(0, vector.lowerBound(4, vector.size, 3))
    }

}
