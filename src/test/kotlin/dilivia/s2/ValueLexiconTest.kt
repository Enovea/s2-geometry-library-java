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
