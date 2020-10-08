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

