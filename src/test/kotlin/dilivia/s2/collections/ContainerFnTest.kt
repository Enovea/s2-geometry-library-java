package dilivia.s2.collections

import junit.framework.TestCase

class IntContainer(val values: List<Int>): Container<Int> {
    override fun size(): Int = values.size
    override fun get(index: Int): Int = values[index]
}

class ContainerFnTest : TestCase() {


    fun testListLowerBound() {
        val vector = IntContainer(listOf(1,2,3,4,5,6))
        assertEquals(2, vector.lowerBound(1, 4, 3))
        assertEquals(4, vector.lowerBound(1, 4, 7))
        assertEquals(4, vector.lowerBound(1, 4, 8))
        assertEquals(2, vector.lowerBound(1, 6, 3))
        assertEquals(6, vector.lowerBound(1, 6, 8))
        assertEquals(4, vector.lowerBound(4, 6, 3))
        assertEquals(6, vector.lowerBound(6 + 1, 6 + 3, 3))
    }

    fun testListLowerBoundEmpty() {
        val vector = listOf<Int>()
        assertEquals(0, vector.lowerBound(1, 4, 3))
        assertEquals(0, vector.lowerBound(1, 4, 7))
        assertEquals(0, vector.lowerBound(1, 4, 8))
        assertEquals(0, vector.lowerBound(1, 6, 3))
        assertEquals(0, vector.lowerBound(1, 6, 8))
        assertEquals(0, vector.lowerBound(4, 6, 3))
    }

}
