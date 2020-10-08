package dilivia.s2.builder

import dilivia.s2.S2GeometryTestCase

class IdSetLexiconTest : S2GeometryTestCase() {

  fun expectIdSet(expected: List<Int>, actual: IdSetLexicon.IdSet) {
    assertEquals(expected.size, actual.size())
    assertEquals(expected, actual.values)
  }

  fun testEmptySet() {
    val lexicon = IdSetLexicon()
    expectIdSet(listOf(), lexicon.idSet(lexicon.add(listOf())))
  }

  fun testSingletonSets() {
    val lexicon = IdSetLexicon()
    assertEquals(5, lexicon.add(listOf( 5 )))
    assertEquals(0, lexicon.add(listOf( 0 )))
    assertEquals(1, lexicon.addSingleton(1))
    val m = Int.MAX_VALUE
    assertEquals(m, lexicon.add(listOf(m)))

    expectIdSet(listOf(0), lexicon.idSet(0))
    expectIdSet(listOf(1), lexicon.idSet(1))
    expectIdSet(listOf(5), lexicon.idSet(5))
    expectIdSet(listOf(m), lexicon.idSet(m))
  }

  fun testSetsAreSorted() {
    val lexicon = IdSetLexicon()
    assertEquals(0.inv(), lexicon.add(listOf(2, 5)))
    assertEquals(1.inv(), lexicon.add(listOf(3, 2, 5)))
    assertEquals(0.inv(), lexicon.add(listOf(5, 2)))
    assertEquals(1.inv(), lexicon.add(listOf(5, 3, 2, 5)))

    expectIdSet(listOf(2, 5), lexicon.idSet(0.inv()))
    expectIdSet(listOf(2, 3, 5), lexicon.idSet(1.inv()))
  }

  fun testClear() {
    val lexicon = IdSetLexicon()
    assertEquals(0.inv(), lexicon.add(listOf(1, 2)))
    assertEquals(1.inv(), lexicon.add(listOf(3, 4)))
    lexicon.clear()
    assertEquals(0.inv(), lexicon.add(listOf(3, 4)))
    assertEquals(1.inv(), lexicon.add(listOf(1, 2)))
  }

}
