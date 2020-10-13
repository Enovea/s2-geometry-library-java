package dilivia.s2.index

import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser
import dilivia.s2.index.S2ShapeIndex.RangeIterator

class S2ShapeIndexRangeIteratorTest : S2GeometryTestCase() {

    fun testNext() {
        // Create an index with one point each on S2CellId faces 0, 1, and 2.
        val index: S2ShapeIndex = S2TextParser.makeIndex("0:0 | 0:90 | 90:0 # #")
        val it = RangeIterator(index)
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
        val empty_it = RangeIterator(empty)
        val non_empty_it = RangeIterator(non_empty)
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
