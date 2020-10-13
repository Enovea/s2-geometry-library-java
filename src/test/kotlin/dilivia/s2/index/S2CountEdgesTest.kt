package dilivia.s2.index

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser

class S2CountEdgesTest : S2GeometryTestCase() {

fun testCountEdgesUpToStopsEarly() {
  val index = S2TextParser.makeIndex(
      "0:0 | 0:1 | 0:2 | 0:3 | 0:4 # 1:0, 1:1 | 1:2, 1:3 | 1:4, 1:5, 1:6 #"
  )
  // Verify the test parameters.
  assertEquals(index.numShapeIds(), 4)
  assertEquals(index.shape(0)?.numEdges, 5)
  assertEquals(index.shape(1)?.numEdges, 1)
  assertEquals(index.shape(2)?.numEdges, 1)
  assertEquals(index.shape(3)?.numEdges, 2)

  assertEquals(S2CountEdges.countEdges(index), 9)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 1), 5)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 5), 5)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 6), 6)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 8), 9)
}

}
